from __future__ import annotations

import argparse
import json
import math
from pathlib import Path


def wilson_interval(successes: int, trials: int, z: float = 1.959963984540054) -> tuple[float, float]:
    if trials < 1 or not 0 <= successes <= trials:
        raise ValueError("invalid binomial counts")
    p = successes / trials
    den = 1.0 + z * z / trials
    centre = (p + z * z / (2.0 * trials)) / den
    half = z * math.sqrt(p * (1.0 - p) / trials + z * z / (4.0 * trials * trials)) / den
    return max(0.0, centre - half), min(1.0, centre + half)


def aggregate(rule: dict, batch_paths: list[Path]) -> dict:
    if rule["status"] != "frozen_before_q5_outcomes":
        raise ValueError("Q5 rule is not frozen")
    expected_n = int(rule["inference"]["replicates_per_cell"])
    alpha = float(rule["inference"]["alpha"])
    cells = {cell["id"]: cell for cell in rule["mandatory_cells"]}
    rows_by_cell: dict[str, dict[int, dict]] = {cell_id: {} for cell_id in cells}

    for path in batch_paths:
        payload = json.loads(path.read_text(encoding="utf-8"))
        if payload.get("schema") != "ttfm_q5_geometry_batch_v0.1":
            raise ValueError(f"unexpected batch schema in {path}")
        cell_id = payload["cell_id"]
        if cell_id not in cells:
            raise ValueError(f"unknown cell in batch: {cell_id}")
        for row in payload["rows"]:
            rep = int(row["replicate"])
            if rep in rows_by_cell[cell_id]:
                raise ValueError(f"duplicate replicate {rep} in {cell_id}")
            rows_by_cell[cell_id][rep] = row

    cell_results = []
    for cell_id, cell in cells.items():
        rep_map = rows_by_cell[cell_id]
        if set(rep_map) != set(range(expected_n)):
            missing = sorted(set(range(expected_n)) - set(rep_map))
            extra = sorted(set(rep_map) - set(range(expected_n)))
            raise ValueError(f"replicate coverage mismatch for {cell_id}: missing={missing[:10]} extra={extra[:10]}")
        rows = [rep_map[i] for i in range(expected_n)]
        m_reject = sum(float(row["ttf_m_p_value"]) <= alpha for row in rows)
        c_reject = sum(float(row["ttf_c_p_value"]) <= alpha for row in rows)
        m_low, m_high = wilson_interval(m_reject, expected_n)
        c_low, c_high = wilson_interval(c_reject, expected_n)
        cell_results.append(
            {
                "id": cell_id,
                "response_mode": cell["response_mode"],
                "geometry_profile": cell["geometry_profile"],
                "roles": list(cell["roles"]),
                "n_replicates": expected_n,
                "TTF_M": {
                    "rejections": int(m_reject),
                    "rejection_rate": float(m_reject / expected_n),
                    "wilson95_low": float(m_low),
                    "wilson95_high": float(m_high),
                    "mean_statistic": float(sum(float(r["ttf_m_statistic"]) for r in rows) / expected_n),
                },
                "TTF_C": {
                    "rejections": int(c_reject),
                    "rejection_rate": float(c_reject / expected_n),
                    "wilson95_low": float(c_low),
                    "wilson95_high": float(c_high),
                    "mean_statistic": float(sum(float(r["ttf_c_statistic"]) for r in rows) / expected_n),
                },
                "geometry_summary": {
                    "mean_train_effective_n": float(sum(float(r["train_effective_n_mean"]) for r in rows) / expected_n),
                    "mean_eval_effective_n": float(sum(float(r["eval_effective_n_mean"]) for r in rows) / expected_n),
                    "minimum_train_effective_n": int(min(int(r["train_effective_n_min"]) for r in rows)),
                    "minimum_eval_effective_n": int(min(int(r["eval_effective_n_min"]) for r in rows)),
                    "mean_train_k": float(sum(float(r["train_k_mean"]) for r in rows) / expected_n),
                    "mean_eval_k": float(sum(float(r["eval_k_mean"]) for r in rows) / expected_n),
                },
            }
        )

    upper = float(rule["inference"]["type1_wilson_upper_ceiling"])
    lower = float(rule["inference"]["power_wilson_lower_floor"])

    def estimand_pass(prefix: str) -> tuple[bool, list[dict]]:
        checks = []
        type_role = "ttf_m_type1" if prefix == "TTF_M" else "ttf_c_type1"
        power_role = "ttf_m_power" if prefix == "TTF_M" else "ttf_c_power"
        for result in cell_results:
            roles = set(result["roles"])
            if type_role in roles:
                observed = float(result[prefix]["wilson95_high"])
                checks.append({
                    "cell_id": result["id"],
                    "role": "type1",
                    "observed_bound": observed,
                    "required": f"<= {upper}",
                    "pass": observed <= upper,
                })
            if power_role in roles:
                observed = float(result[prefix]["wilson95_low"])
                checks.append({
                    "cell_id": result["id"],
                    "role": "power",
                    "observed_bound": observed,
                    "required": f">= {lower}",
                    "pass": observed >= lower,
                })
        return all(check["pass"] for check in checks), checks

    m_pass, m_checks = estimand_pass("TTF_M")
    c_pass, c_checks = estimand_pass("TTF_C")
    return {
        "schema": "ttfm_q5_geometry_stress_result_v0.1",
        "status": "formal_q5_geometry_stress_complete",
        "rule_schema": rule["schema"],
        "TTF_M_Q5_pass": bool(m_pass),
        "TTF_C_Q5_pass": bool(c_pass),
        "joint_Q5_pass": bool(m_pass and c_pass),
        "TTF_M_checks": m_checks,
        "TTF_C_checks": c_checks,
        "cells": cell_results,
        "claim_firewall": rule["claim_firewall"],
        "next_if_pass": "independent_empirical_observation_support_and_application_specific_measurement_contract",
        "next_if_fail": "preserve_Q5_fail_and_diagnose_geometry_limit_without_retuning_this_gate"
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--rule", required=True)
    parser.add_argument("--batch-dir", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()
    rule = json.loads(Path(args.rule).read_text(encoding="utf-8"))
    batch_paths = sorted(Path(args.batch_dir).rglob("*.json"))
    if not batch_paths:
        raise ValueError("no Q5 batch files found")
    result = aggregate(rule, batch_paths)
    out = Path(args.output)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps({
        "TTF_M_Q5_pass": result["TTF_M_Q5_pass"],
        "TTF_C_Q5_pass": result["TTF_C_Q5_pass"],
        "joint_Q5_pass": result["joint_Q5_pass"],
    }, sort_keys=True))


if __name__ == "__main__":
    main()
