from __future__ import annotations

import argparse
import json
from pathlib import Path

from ttf.precision import wilson_interval


def aggregate(protocol: dict, development: dict, q5_1: dict, batch_paths: list[Path]) -> dict:
    if protocol.get("status") != "frozen_before_q5_2_development_aggregate_and_formal_outcomes":
        raise ValueError("Q5.2 conditional formal protocol is not frozen")
    if development.get("schema") != "ttfm_q5_2_geometry_transport_development_result_v0.1":
        raise ValueError("unexpected development result schema")
    candidate = development.get("candidate_selected")
    if candidate is None:
        raise ValueError("formal aggregation forbidden because development selected no candidate")
    allowed = set(protocol["candidate_handoff"]["allowed_values"])
    if candidate not in allowed:
        raise ValueError("selected candidate outside frozen formal protocol")
    checks = development.get("candidate_checks", {})
    if candidate not in checks or checks[candidate].get("development_pass") is not True:
        raise ValueError("selected candidate did not pass frozen development gate")
    if q5_1.get("TTF_M_Q5_1_pass") is not True:
        raise ValueError("inherited TTF-M Q5.1 qualification is not PASS")
    if q5_1.get("TTF_C_Q5_1_pass") is not False:
        raise ValueError("Q5.1 TTF-C failure history drifted")

    expected_n = int(protocol["formal_worlds"]["worlds_per_cell"])
    alpha = float(protocol["formal_worlds"]["alpha"])
    cells = {cell["id"]: cell for cell in protocol["formal_cells"]}
    rows_by_cell: dict[str, dict[int, dict]] = {cell_id: {} for cell_id in cells}

    for path in batch_paths:
        payload = json.loads(path.read_text(encoding="utf-8"))
        if payload.get("schema") != "ttfm_q5_2_formal_batch_v0.1":
            raise ValueError(f"unexpected formal batch schema in {path}")
        if payload.get("candidate") != candidate:
            raise ValueError("formal batch candidate drift")
        cell_id = payload.get("cell_id")
        if cell_id not in cells:
            raise ValueError(f"unknown formal cell: {cell_id}")
        for row in payload["rows"]:
            if row.get("candidate") != candidate:
                raise ValueError("row candidate drift")
            rep = int(row["replicate"])
            if rep in rows_by_cell[cell_id]:
                raise ValueError(f"duplicate replicate {rep} in {cell_id}")
            rows_by_cell[cell_id][rep] = row

    gate = protocol["precision_gate"]
    cell_results = []
    all_pass = True
    for cell_id, cell in cells.items():
        rep_map = rows_by_cell[cell_id]
        expected = set(range(expected_n))
        if set(rep_map) != expected:
            missing = sorted(expected - set(rep_map))
            extra = sorted(set(rep_map) - expected)
            raise ValueError(f"formal replicate coverage mismatch for {cell_id}: missing={missing[:10]} extra={extra[:10]}")
        rows = [rep_map[i] for i in range(expected_n)]
        rejects = sum(float(row["p_value"]) <= alpha for row in rows)
        interval = wilson_interval(rejects, expected_n)
        role = str(cell["role"])
        if role == "type1":
            passed = interval.high <= float(gate["type1_wilson_upper_ceiling"])
            boundary = "upper95"
            required = f"<= {float(gate['type1_wilson_upper_ceiling'])}"
        elif role == "power":
            passed = interval.low >= float(gate["power_wilson_lower_floor"])
            boundary = "lower95"
            required = f">= {float(gate['power_wilson_lower_floor'])}"
        else:
            raise ValueError(f"unknown formal role: {role}")
        all_pass = bool(all_pass and passed)
        cell_results.append(
            {
                "id": cell_id,
                "response_mode": cell["response_mode"],
                "geometry_profile": cell["geometry_profile"],
                "role": role,
                "n_replicates": expected_n,
                "rejections": int(rejects),
                "rejection_rate": float(rejects / expected_n),
                "wilson95": interval.to_dict(),
                "decision_boundary": boundary,
                "required": required,
                "pass": bool(passed),
                "mean_statistic": float(sum(float(row["statistic"]) for row in rows) / expected_n),
                "mean_train_effective_edge_count": float(sum(float(row["mean_train_effective_edge_count"]) for row in rows) / expected_n),
                "minimum_train_effective_edge_count": float(min(float(row["min_train_effective_edge_count"]) for row in rows)),
            }
        )

    ttf_m_pass = bool(q5_1["TTF_M_Q5_1_pass"])
    ttf_c_pass = bool(all_pass)
    joint_pass = bool(ttf_m_pass and ttf_c_pass)
    return {
        "schema": "ttfm_q5_2_formal_result_v0.1",
        "status": "formal_complete",
        "candidate": candidate,
        "TTF_M_Q5_2_pass": ttf_m_pass,
        "TTF_C_Q5_2_pass": ttf_c_pass,
        "joint_Q5_2_pass": joint_pass,
        "ttf_m_source": "inherited immutable Q5.1 PASS; TTF-M estimator is unchanged by Q5.2",
        "cells": cell_results,
        "precision_gate": gate,
        "q5_remains_immutable_fail": True,
        "q5_1_remains_immutable_fail": True,
        "formal_candidate_may_not_be_retuned": True,
        "empirical_opening_permitted": bool(joint_pass),
        "empirical_next_gate": (
            protocol["terminal_rules"]["empirical_next_gate_after_pass"]
            if joint_pass
            else None
        ),
        "claim_firewall": protocol["claim_firewall"],
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--protocol", required=True)
    parser.add_argument("--development-result", required=True)
    parser.add_argument("--q5-1-result", required=True)
    parser.add_argument("--batch-dir", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    protocol = json.loads(Path(args.protocol).read_text(encoding="utf-8"))
    development = json.loads(Path(args.development_result).read_text(encoding="utf-8"))
    q5_1 = json.loads(Path(args.q5_1_result).read_text(encoding="utf-8"))
    batch_paths = sorted(Path(args.batch_dir).rglob("*.json"))
    if not batch_paths:
        raise ValueError("no Q5.2 formal batches found")
    result = aggregate(protocol, development, q5_1, batch_paths)
    out = Path(args.output)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps({
        "candidate": result["candidate"],
        "TTF_C_Q5_2_pass": result["TTF_C_Q5_2_pass"],
        "joint_Q5_2_pass": result["joint_Q5_2_pass"],
    }, sort_keys=True))


if __name__ == "__main__":
    main()
