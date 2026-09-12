from __future__ import annotations

import argparse
import json
from pathlib import Path


AUDIT_BASELINE = "q5_1_uniform"


def aggregate(rule: dict, batch_paths: list[Path]) -> dict:
    if rule["status"] != "frozen_before_q5_2_transport_development_outcomes":
        raise ValueError("Q5.2 transport development rule is not frozen")
    expected_n = int(rule["development_inference"]["worlds_per_cell"])
    alpha = float(rule["development_inference"]["alpha"])
    cells = {cell["id"]: cell for cell in rule["development_cells"]}
    modes = [AUDIT_BASELINE] + list(rule["candidate_preference_order"])
    rows_by_cell: dict[str, dict[int, dict]] = {cell_id: {} for cell_id in cells}

    for path in batch_paths:
        payload = json.loads(path.read_text(encoding="utf-8"))
        if payload.get("schema") != "ttfm_q5_2_geometry_transport_development_batch_v0.1":
            raise ValueError(f"unexpected batch schema in {path}")
        cell_id = payload["cell_id"]
        if cell_id not in cells:
            raise ValueError(f"unknown cell: {cell_id}")
        for row in payload["rows"]:
            rep = int(row["replicate"])
            if rep in rows_by_cell[cell_id]:
                raise ValueError(f"duplicate replicate {rep} in {cell_id}")
            if set(row["modes"]) != set(modes):
                raise ValueError(f"mode drift in {cell_id} replicate {rep}")
            rows_by_cell[cell_id][rep] = row

    cell_results = []
    for cell_id, cell in cells.items():
        rep_map = rows_by_cell[cell_id]
        if set(rep_map) != set(range(expected_n)):
            missing = sorted(set(range(expected_n)) - set(rep_map))
            extra = sorted(set(rep_map) - set(range(expected_n)))
            raise ValueError(
                f"replicate coverage mismatch for {cell_id}: missing={missing[:10]} extra={extra[:10]}"
            )
        rows = [rep_map[i] for i in range(expected_n)]
        mode_summary = {}
        for mode in modes:
            rejects = sum(float(row["modes"][mode]["p_value"]) <= alpha for row in rows)
            mean_stat = sum(float(row["modes"][mode]["statistic"]) for row in rows) / expected_n
            mean_eff = sum(
                float(row["modes"][mode]["mean_train_effective_edge_count"])
                for row in rows
            ) / expected_n
            min_eff = min(
                float(row["modes"][mode]["min_train_effective_edge_count"])
                for row in rows
            )
            mode_summary[mode] = {
                "rejections": int(rejects),
                "rejection_rate": float(rejects / expected_n),
                "mean_statistic": float(mean_stat),
                "mean_train_effective_edge_count": float(mean_eff),
                "minimum_train_effective_edge_count": float(min_eff),
            }
        baseline_rate = float(mode_summary[AUDIT_BASELINE]["rejection_rate"])
        for mode in rule["candidate_preference_order"]:
            mode_summary[mode]["rejection_rate_delta_from_q5_1"] = float(
                mode_summary[mode]["rejection_rate"] - baseline_rate
            )
        cell_results.append(
            {
                "id": cell_id,
                "response_mode": cell["response_mode"],
                "geometry_profile": cell["geometry_profile"],
                "role": cell["role"],
                "n_replicates": expected_n,
                "modes": mode_summary,
            }
        )

    thresholds = rule["selection_rule"]
    candidate_checks = {}
    for mode in rule["candidate_preference_order"]:
        checks = []
        for result in cell_results:
            rate = float(result["modes"][mode]["rejection_rate"])
            role = result["role"]
            if role == "type1":
                required = f"<= {float(thresholds['type1_max_rate'])}"
                passed = rate <= float(thresholds["type1_max_rate"])
            elif role == "shared_mismatch_matched_power":
                required = f">= {float(thresholds['shared_mismatch_matched_min_power'])}"
                passed = rate >= float(thresholds["shared_mismatch_matched_min_power"])
            elif role == "shared_mismatch_shifted_power":
                required = f">= {float(thresholds['shared_mismatch_shifted_min_power'])}"
                passed = rate >= float(thresholds["shared_mismatch_shifted_min_power"])
            elif role == "shared_relation_matched_power":
                required = f">= {float(thresholds['shared_relation_matched_min_power'])}"
                passed = rate >= float(thresholds["shared_relation_matched_min_power"])
            elif role == "shared_relation_shifted_power":
                required = f">= {float(thresholds['shared_relation_shifted_min_power'])}"
                passed = rate >= float(thresholds["shared_relation_shifted_min_power"])
            else:
                raise ValueError(f"unknown development role: {role}")
            checks.append(
                {
                    "cell_id": result["id"],
                    "role": role,
                    "observed_rate": rate,
                    "required": required,
                    "pass": bool(passed),
                }
            )
        candidate_checks[mode] = {
            "development_pass": bool(all(check["pass"] for check in checks)),
            "checks": checks,
        }

    selected = None
    for mode in rule["candidate_preference_order"]:
        if candidate_checks[mode]["development_pass"]:
            selected = mode
            break

    return {
        "schema": "ttfm_q5_2_geometry_transport_development_result_v0.1",
        "status": "development_complete_nonqualifying",
        "formal_q5_2_qualified": False,
        "candidate_selected": selected,
        "selection_rule": rule["selection_rule"],
        "candidate_checks": candidate_checks,
        "cells": cell_results,
        "baseline_mode": AUDIT_BASELINE,
        "q5_remains_immutable_fail": True,
        "q5_1_remains_immutable_fail": True,
        "next_step": (
            rule["next_if_candidate_selected"]
            if selected is not None
            else rule["next_if_no_candidate"]
        ),
        "claim_firewall": rule["claim_firewall"],
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
        raise ValueError("no Q5.2 transport development batches found")
    result = aggregate(rule, batch_paths)
    out = Path(args.output)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(
        json.dumps(
            {
                "candidate_selected": result["candidate_selected"],
                "formal_q5_2_qualified": result["formal_q5_2_qualified"],
            },
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()
