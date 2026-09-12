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


def _rate_block(successes: int, trials: int) -> dict:
    low, high = wilson_interval(successes, trials)
    return {
        "count": int(successes),
        "rate": float(successes / trials),
        "wilson95_low": float(low),
        "wilson95_high": float(high),
    }


def aggregate(rule: dict, batch_paths: list[Path]) -> dict:
    if rule["status"] != "frozen_before_q4_1_outcomes":
        raise ValueError("Q4.1 rule is not frozen")
    expected_n = int(rule["inference"]["replicates_per_cell"])
    upper = float(rule["inference"]["type1_wilson_upper_max"])
    lower = float(rule["inference"]["power_wilson_lower_min"])
    cells = {cell["id"]: cell for cell in rule["mandatory_cells"]}
    rows_by_cell: dict[str, dict[int, dict]] = {cell_id: {} for cell_id in cells}

    for path in batch_paths:
        payload = json.loads(path.read_text(encoding="utf-8"))
        if payload.get("schema") != "ttfm_q4_1_strict_private_batch_v0.1":
            raise ValueError(f"unexpected Q4.1 batch schema in {path}")
        cell_id = payload["cell_id"]
        if cell_id not in cells:
            raise ValueError(f"unknown Q4.1 cell: {cell_id}")
        for row in payload["rows"]:
            rep = int(row["replicate"])
            if rep in rows_by_cell[cell_id]:
                raise ValueError(f"duplicate Q4.1 replicate {rep} in {cell_id}")
            rows_by_cell[cell_id][rep] = row

    cell_results = []
    checks = []
    metric_for_role = {
        "relation": "relation_detection",
        "breakdown": "breakdown_label",
        "recoupling": "recoupling_label",
    }

    for cell_id, cell in cells.items():
        rep_map = rows_by_cell[cell_id]
        if set(rep_map) != set(range(expected_n)):
            missing = sorted(set(range(expected_n)) - set(rep_map))
            extra = sorted(set(rep_map) - set(range(expected_n)))
            raise ValueError(
                f"Q4.1 replicate coverage mismatch {cell_id}: "
                f"missing={missing[:10]} extra={extra[:10]}"
            )
        rows = [rep_map[i] for i in range(expected_n)]
        relation = sum(bool(r["relation_detected"]) for r in rows)
        breakdown = sum(bool(r["breakdown_label"]) for r in rows)
        recoupling = sum(bool(r["recoupling_label"]) for r in rows)
        result = {
            "id": cell_id,
            "mode": cell["mode"],
            "roles": list(cell["roles"]),
            "noise_sd": float(cell["noise_sd"]),
            "n_replicates": int(expected_n),
            "relation_detection": _rate_block(relation, expected_n),
            "breakdown_label": _rate_block(breakdown, expected_n),
            "recoupling_label": _rate_block(recoupling, expected_n),
            "mean_front_center": float(sum(float(r["front_center"]) for r in rows) / expected_n),
            "mean_oriented_delta": float(sum(float(r["mean_oriented_delta"]) for r in rows) / expected_n),
            "mean_positive_fraction": float(sum(float(r["positive_fraction"]) for r in rows) / expected_n),
            "mean_negative_fraction": float(sum(float(r["negative_fraction"]) for r in rows) / expected_n),
            "mean_ttf_c_statistic": float(sum(float(r["ttf_c_statistic"]) for r in rows) / expected_n),
        }
        cell_results.append(result)

        for role in cell["roles"]:
            if role == "characterization_only":
                continue
            prefix, kind = role.rsplit("_", 1)
            if prefix not in metric_for_role:
                raise ValueError(f"unknown Q4.1 role metric: {role}")
            block = result[metric_for_role[prefix]]
            if kind == "type1":
                observed = float(block["wilson95_high"])
                passed = observed <= upper
                required = f"<= {upper}"
            elif kind == "power":
                observed = float(block["wilson95_low"])
                passed = observed >= lower
                required = f">= {lower}"
            else:
                raise ValueError(f"unknown Q4.1 role kind: {role}")
            checks.append(
                {
                    "cell_id": cell_id,
                    "role": role,
                    "observed_bound": observed,
                    "required": required,
                    "pass": bool(passed),
                }
            )

    q4_1_pass = bool(all(check["pass"] for check in checks))
    return {
        "schema": "ttfm_q4_1_strict_private_result_v0.1",
        "status": "formal_q4_1_strict_private_complete",
        "rule_schema": rule["schema"],
        "Q4_1_pass": q4_1_pass,
        "checks": checks,
        "cells": cell_results,
        "characterization_cells_excluded_from_gate": [
            x["id"] for x in cell_results if "characterization_only" in x["roles"]
        ],
        "history_firewall": rule["history_firewall"],
        "claim_firewall": rule["claim_firewall"],
        "next_required_if_pass": "Q5_heterogeneous_geometry_stress",
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
        raise ValueError("no Q4.1 batch files found")
    result = aggregate(rule, batch_paths)
    out = Path(args.output)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps({"Q4_1_pass": result["Q4_1_pass"]}, sort_keys=True))


if __name__ == "__main__":
    main()
