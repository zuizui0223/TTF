#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
from pathlib import Path

from ttf.genetic_geometry_io import sha256_path


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--input-dir", type=Path, required=True)
    ap.add_argument("--execution-plan", type=Path, required=True)
    ap.add_argument("--phase3-rule", type=Path, required=True)
    ap.add_argument("--retention", type=float, required=True)
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()

    plan = json.loads(args.execution_plan.read_text())
    rule = json.loads(args.phase3_rule.read_text())
    if plan.get("schema") != "ttf_genetic_phylogatr_gate_d_fragility_execution_plan_v0.1":
        raise RuntimeError("fragility execution plan schema drift")
    if rule.get("schema") != "ttf_genetic_phylogatr_phase3_gate_d_rule_v0.1":
        raise RuntimeError("formal Phase-3 rule schema drift")
    if sha256_path(args.phase3_rule) != plan["phase3_rule_sha256"]:
        raise RuntimeError("formal Phase-3 rule provenance drift")

    retention = float(args.retention)
    levels = [level for level in plan["levels"] if float(level["retention_fraction"]) == retention]
    if len(levels) != 1 or levels[0]["status"] != "SYNTHETIC_DIAGNOSTIC_RUNNABLE":
        raise RuntimeError("requested fragility level is not runnable")
    level = levels[0]
    expected_labels = tuple(rule["synthetic_world"]["private_reference_configurations"].keys())
    expected_n = int(rule["qualification"]["reference_worlds_per_configuration"])
    grouped: dict[str, dict[int, tuple[float, float]]] = {label: {} for label in expected_labels}

    for path in sorted(args.input_dir.rglob("*.json")):
        shard = json.loads(path.read_text())
        if shard.get("schema") != "ttf_genetic_phylogatr_gate_d_fragility_reference_shard_v0.1":
            continue
        if float(shard["retention_fraction"]) != retention:
            continue
        if shard.get("geometry_fingerprint_sha256") != level["geometry_fingerprint_sha256"]:
            raise RuntimeError(f"fragility reference geometry fingerprint drift in {path}")
        for key in (
            "confirmatory_sequence_identity_opened",
            "confirmatory_pairwise_genetic_distances_opened",
            "confirmatory_ttf_statistic_opened",
            "formal_gate_d_decision_authority",
            "phase4_identity_opening_authority",
        ):
            if shard.get(key) is not False:
                raise RuntimeError(f"fragility reference authority/firewall drift in {path}: {key}")
        label = str(shard["configuration"])
        if label not in grouped:
            raise RuntimeError(f"unexpected fragility reference configuration {label}")
        for row in shard["rows"]:
            replicate = int(row["replicate"])
            if replicate in grouped[label]:
                raise RuntimeError(f"duplicate fragility reference replicate {label}:{replicate}")
            grouped[label][replicate] = (
                float(row["training_strength"]),
                float(row["statistic"]),
            )

    references: dict[str, dict[str, list[float]]] = {}
    for label in expected_labels:
        rows = grouped[label]
        if sorted(rows) != list(range(expected_n)):
            missing = sorted(set(range(expected_n)) - set(rows))[:10]
            raise RuntimeError(
                f"incomplete fragility reference configuration {label}; first missing={missing}"
            )
        references[label] = {
            "training_strength": [rows[index][0] for index in range(expected_n)],
            "statistic": [rows[index][1] for index in range(expected_n)],
        }

    qualification = rule["qualification"]
    out = {
        "schema": "ttf_genetic_phylogatr_gate_d_fragility_references_v0.1",
        "status": "DIAGNOSTIC_PRIVATE_REFERENCE_FAMILIES_COMPLETE",
        "retention_fraction": retention,
        "retention_label": level["retention_label"],
        "geometry_fingerprint_sha256": level["geometry_fingerprint_sha256"],
        "reference_worlds_per_configuration": expected_n,
        "profile_draws": int(qualification["profile_strength_draws"]),
        "calibration_statistic_draws": int(qualification["calibration_statistic_draws"]),
        "configuration_order": list(expected_labels),
        "references": references,
        "confirmatory_sequence_identity_opened": False,
        "confirmatory_pairwise_genetic_distances_opened": False,
        "confirmatory_ttf_statistic_opened": False,
        "formal_gate_d_decision_authority": False,
        "phase4_identity_opening_authority": False,
    }
    if out["profile_draws"] + out["calibration_statistic_draws"] != expected_n:
        raise RuntimeError("fragility profile/calibration split does not cover references exactly")
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(out, indent=2, sort_keys=True) + "\n")
    print(json.dumps({key: value for key, value in out.items() if key != "references"}, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
