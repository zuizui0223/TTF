#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
from pathlib import Path

from ttf.precision import wilson_interval


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--input-dir", type=Path, required=True)
    ap.add_argument("--self-rule", type=Path, required=True)
    ap.add_argument("--phase3-authorization", type=Path, required=True)
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()

    rule = json.loads(args.self_rule.read_text())
    auth = json.loads(args.phase3_authorization.read_text())
    if rule.get("schema") != "ttf_genetic_phylogatr_phase3_self_detectability_rule_v0.1":
        raise RuntimeError("fresh self rule schema drift")
    if auth.get("schema") != "ttf_genetic_phylogatr_phase3_gate_d_authorization_v0.1":
        raise RuntimeError("fresh Phase-3 authorization schema drift")
    if any(value is not False for value in rule["outcome_firewall"].values()):
        raise RuntimeError("fresh self-detectability outcome firewall is open")

    syn = rule["synthetic_worlds"]
    expected = {
        "null": int(syn["null_evaluation_worlds"]),
        "private_A2": int(syn["private_A2_evaluation_worlds"]),
    }
    grouped: dict[str, dict[int, float]] = {cell: {} for cell in expected}
    for path in sorted(args.input_dir.rglob("*.json")):
        payload = json.loads(path.read_text())
        if payload.get("schema") != "ttf_genetic_phylogatr_phase3_self_evaluation_shard_v0.1":
            continue
        if payload.get("geometry_fingerprint_sha256") != auth["geometry_fingerprint_sha256"]:
            raise RuntimeError(f"fresh self-evaluation geometry drift in {path}")
        if payload.get("confirmatory_sequence_identity_opened") is not False:
            raise RuntimeError(f"fresh identity firewall drift in {path}")
        if payload.get("confirmatory_pairwise_genetic_distances_opened") is not False:
            raise RuntimeError(f"fresh distance firewall drift in {path}")
        cell = str(payload["cell"])
        if cell not in grouped:
            raise RuntimeError(f"unexpected fresh self cell {cell}")
        for row in payload["rows"]:
            replicate = int(row["replicate"])
            if replicate in grouped[cell]:
                raise RuntimeError(f"duplicate fresh self replicate {cell}:{replicate}")
            grouped[cell][replicate] = float(row["p_value"])

    alpha = float(rule["inference"]["alpha"])
    cells: dict[str, dict[str, float | int]] = {}
    for cell, n_worlds in expected.items():
        if sorted(grouped[cell]) != list(range(n_worlds)):
            missing = sorted(set(range(n_worlds)) - set(grouped[cell]))[:10]
            raise RuntimeError(f"incomplete fresh self cell {cell}; first missing={missing}")
        rejections = sum(grouped[cell][index] <= alpha for index in range(n_worlds))
        interval = wilson_interval(rejections, n_worlds)
        cells[cell] = {
            "n_worlds": n_worlds,
            "rejections": int(rejections),
            "rejection_rate": float(rejections / n_worlds),
            "wilson95_low": float(interval.low),
            "wilson95_high": float(interval.high),
        }

    qualification = rule["qualification"]
    type1_ceiling = float(qualification["type1_wilson95_upper_ceiling"])
    power_floor = float(qualification["private_A2_wilson95_lower_floor"])
    type1_pass = float(cells["null"]["wilson95_high"]) <= type1_ceiling
    power_pass = float(cells["private_A2"]["wilson95_low"]) >= power_floor
    passed = bool(type1_pass and power_pass)

    out = {
        "schema": "ttf_genetic_phylogatr_phase3_self_qualification_v0.1",
        "status": "PASS" if passed else "SELF_DETECTABILITY_NOT_QUALIFIED",
        "geometry_fingerprint_sha256": auth["geometry_fingerprint_sha256"],
        "cells": cells,
        "type1_gate": {
            "wilson95_upper_ceiling": type1_ceiling,
            "observed_wilson95_upper": float(cells["null"]["wilson95_high"]),
            "pass": bool(type1_pass),
        },
        "power_gate": {
            "wilson95_lower_floor": power_floor,
            "observed_wilson95_lower": float(cells["private_A2"]["wilson95_low"]),
            "pass": bool(power_pass),
        },
        "passed": passed,
        "confirmatory_sequence_identity_opened": False,
        "confirmatory_pairwise_genetic_distances_opened": False,
        "confirmatory_ttf_statistic_opened": False,
        "interpretation_contract": {
            "if_future_cross_species_primary_positive": "transferable_place_component; this self gate is not required for the positive decision",
            "if_future_cross_species_primary_negative_and_self_pass": "lineage_conditioned_spatial_structure_within_tested_domain",
            "if_future_cross_species_primary_negative_and_self_not_pass": "NOT_EVALUABLE_FOR_LINEAGE_CONDITIONING",
        },
        "claim_boundary": (
            "This receipt qualifies only the within-species self-detectability diagnostic on the exact fresh Phase-2 geometry. "
            "It does not open nucleotide identity and does not constitute an empirical genetic result."
        ),
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(out, indent=2, sort_keys=True) + "\n")
    print(json.dumps(out, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
