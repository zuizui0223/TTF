#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
from pathlib import Path


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
    if auth.get("status") != "authorize_frozen_fresh_phylogatr_phase3_gate_d":
        raise RuntimeError("fresh phase-3 Gate-D is not authorized")

    expected_labels = tuple(rule["synthetic_world"]["private_reference_configurations"].keys())
    expected_n = int(rule["qualification"]["reference_worlds_per_configuration"])
    grouped: dict[str, dict[int, tuple[float, float]]] = {label: {} for label in expected_labels}

    files = sorted(args.input_dir.rglob("*.json"))
    if not files:
        raise RuntimeError("no fresh phase-3 reference shard files found")
    for path in files:
        shard = json.loads(path.read_text())
        if shard.get("schema") != "ttf_genetic_phylogatr_phase3_reference_shard_v0.1":
            continue
        if shard.get("confirmatory_sequence_identity_opened") is not False:
            raise RuntimeError(f"fresh nucleotide-identity firewall drift in {path}")
        if shard.get("confirmatory_pairwise_genetic_distances_opened") is not False:
            raise RuntimeError(f"fresh distance firewall drift in {path}")
        if shard.get("geometry_fingerprint_sha256") != auth["geometry_fingerprint_sha256"]:
            raise RuntimeError(f"fresh reference geometry fingerprint drift in {path}")
        label = str(shard["configuration"])
        if label not in grouped:
            raise RuntimeError(f"unexpected fresh reference configuration {label}")
        for row in shard["rows"]:
            replicate = int(row["replicate"])
            if replicate in grouped[label]:
                raise RuntimeError(f"duplicate fresh reference replicate {label}:{replicate}")
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
                f"incomplete fresh reference configuration {label}; first missing={missing}"
            )
        references[label] = {
            "training_strength": [rows[i][0] for i in range(expected_n)],
            "statistic": [rows[i][1] for i in range(expected_n)],
        }

    qualification = rule["qualification"]
    out = {
        "schema": "ttf_genetic_phylogatr_phase3_references_v0.1",
        "status": "ordered_private_reference_families_complete",
        "geometry_fingerprint_sha256": auth["geometry_fingerprint_sha256"],
        "reference_worlds_per_configuration": expected_n,
        "profile_draws": int(qualification["profile_strength_draws"]),
        "calibration_statistic_draws": int(qualification["calibration_statistic_draws"]),
        "configuration_order": list(expected_labels),
        "references": references,
        "confirmatory_sequence_identity_opened": False,
        "confirmatory_pairwise_genetic_distances_opened": False,
        "qualification_claim_made": False,
    }
    if out["profile_draws"] + out["calibration_statistic_draws"] != expected_n:
        raise RuntimeError("fresh profile/calibration split does not cover references exactly")
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(out, indent=2, sort_keys=True) + "\n")
    print(json.dumps({key: value for key, value in out.items() if key != "references"}, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
