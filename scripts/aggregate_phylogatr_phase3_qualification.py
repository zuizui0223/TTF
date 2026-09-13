#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
from collections import Counter
from pathlib import Path

from ttf.precision import wilson_interval


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--input-dir", type=Path, required=True)
    ap.add_argument("--phase3-rule", type=Path, required=True)
    ap.add_argument("--authorization", type=Path, required=True)
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()

    rule = json.loads(args.phase3_rule.read_text())
    auth = json.loads(args.authorization.read_text())
    if rule.get("schema") != "ttf_genetic_phylogatr_phase3_gate_d_rule_v0.1":
        raise RuntimeError("fresh phase-3 rule schema drift")
    if auth.get("schema") != "ttf_genetic_phylogatr_phase3_gate_d_authorization_v0.1":
        raise RuntimeError("fresh phase-3 authorization schema drift")

    qualification = rule["qualification"]
    expected_cells = [
        (float(cell[0]), float(cell[1]))
        for cell in qualification["mandatory_primary_cells"]
    ]
    expected_n = int(qualification["observed_worlds_per_cell"])
    alpha = float(qualification["alpha"])
    grouped: dict[tuple[float, float], dict[int, dict]] = {cell: {} for cell in expected_cells}

    files = sorted(args.input_dir.rglob("*.json"))
    if not files:
        raise RuntimeError("no fresh phase-3 observed shard files found")
    for path in files:
        shard = json.loads(path.read_text())
        if shard.get("schema") != "ttf_genetic_phylogatr_phase3_observed_shard_v0.1":
            continue
        if shard.get("confirmatory_sequence_identity_opened") is not False:
            raise RuntimeError(f"fresh identity firewall drift in {path}")
        if shard.get("confirmatory_pairwise_genetic_distances_opened") is not False:
            raise RuntimeError(f"fresh distance firewall drift in {path}")
        if shard.get("geometry_fingerprint_sha256") != auth["geometry_fingerprint_sha256"]:
            raise RuntimeError(f"fresh observed geometry fingerprint drift in {path}")
        key = (float(shard["shared_fraction"]), float(shard["residual_amplitude"]))
        if key not in grouped:
            raise RuntimeError(f"unexpected fresh observed cell {key}")
        for row in shard["rows"]:
            replicate = int(row["replicate"])
            if replicate in grouped[key]:
                raise RuntimeError(f"duplicate fresh observed replicate {key}:{replicate}")
            grouped[key][replicate] = row

    cells: list[dict] = []
    for shared, amplitude in expected_cells:
        rows = grouped[(shared, amplitude)]
        if sorted(rows) != list(range(expected_n)):
            missing = sorted(set(range(expected_n)) - set(rows))[:10]
            raise RuntimeError(f"incomplete fresh observed cell {(shared, amplitude)}; first missing={missing}")
        ordered = [rows[i] for i in range(expected_n)]
        rejections = sum(float(row["p_value"]) <= alpha for row in ordered)
        interval = wilson_interval(rejections, expected_n)
        selected_counts = Counter(
            "+".join(row["selected_configurations"]) for row in ordered
        )
        cells.append(
            {
                "shared_fraction": shared,
                "residual_amplitude": amplitude,
                "n_worlds": expected_n,
                "rejections": int(rejections),
                "rejection_rate": float(rejections / expected_n),
                "wilson95_low": float(interval.low),
                "wilson95_high": float(interval.high),
                "selected_pair_counts": dict(sorted(selected_counts.items())),
            }
        )

    type1_cells = [cell for cell in cells if cell["shared_fraction"] == 0.0]
    shared_a2 = next(
        cell
        for cell in cells
        if cell["shared_fraction"] == 1.0 and cell["residual_amplitude"] == 2.0
    )
    type1_ceiling = float(qualification["type1_wilson95_upper_ceiling"])
    power_floor = float(qualification["shared_A2_wilson95_lower_floor"])
    max_type1_upper = max(float(cell["wilson95_high"]) for cell in type1_cells)
    type1_pass = all(float(cell["wilson95_high"]) <= type1_ceiling for cell in type1_cells)
    power_pass = float(shared_a2["wilson95_low"]) >= power_floor
    passed = bool(type1_pass and power_pass)

    out = {
        "schema": "ttf_genetic_phylogatr_phase3_qualification_v0.1",
        "status": "PASS" if passed else "NOT_EVALUABLE",
        "geometry_fingerprint_sha256": auth["geometry_fingerprint_sha256"],
        "cells": cells,
        "type1_gate": {
            "wilson95_upper_ceiling": type1_ceiling,
            "max_observed_wilson95_upper": max_type1_upper,
            "pass": bool(type1_pass),
        },
        "power_gate": {
            "cell": {"shared_fraction": 1.0, "residual_amplitude": 2.0},
            "wilson95_lower_floor": power_floor,
            "observed_wilson95_lower": float(shared_a2["wilson95_low"]),
            "pass": bool(power_pass),
        },
        "passed": passed,
        "failure_interpretation": qualification["failure_interpretation"],
        "phase4_identity_opening_eligible": passed,
        "confirmatory_sequence_identity_opened": False,
        "confirmatory_pairwise_genetic_distances_opened": False,
        "confirmatory_ttf_statistic_opened": False,
        "claim_boundary": (
            "PASS authorizes only a separate explicit Phase-4 nucleotide-identity opening receipt "
            "for this exact fresh dataset and split. PASS is not itself an empirical biological result. "
            "NOT_EVALUABLE is not evidence against transferable phylogeography."
        ),
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(out, indent=2, sort_keys=True) + "\n")
    print(json.dumps(out, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
