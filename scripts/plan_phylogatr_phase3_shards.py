#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
from pathlib import Path


def _ranges(total: int, size: int) -> list[tuple[int, int]]:
    if total < 1 or size < 1:
        raise ValueError("total and shard size must be positive")
    return [(start, min(start + size, total)) for start in range(0, total, size)]


def main() -> int:
    ap = argparse.ArgumentParser()
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
    if any(value is not False for value in auth["outcome_firewall"].values()):
        raise RuntimeError("fresh phase-3 empirical outcome firewall is open")

    qualification = rule["qualification"]
    reference_n = int(qualification["reference_worlds_per_configuration"])
    observed_n = int(qualification["observed_worlds_per_cell"])
    reference_ranges = _ranges(reference_n, int(auth["execution"]["reference_shard_size"]))
    observed_ranges = _ranges(observed_n, int(auth["execution"]["observed_shard_size"]))

    reference_jobs = [
        {"configuration": label, "start": start, "stop": stop}
        for label in rule["synthetic_world"]["private_reference_configurations"]
        for start, stop in reference_ranges
    ]
    observed_jobs = [
        {
            "shared_fraction": float(cell[0]),
            "residual_amplitude": float(cell[1]),
            "start": start,
            "stop": stop,
        }
        for cell in qualification["mandatory_primary_cells"]
        for start, stop in observed_ranges
    ]

    if len(reference_jobs) != int(auth["execution"]["reference_expected_jobs"]):
        raise RuntimeError("fresh reference job-count drift")
    if len(observed_jobs) != int(auth["execution"]["observed_expected_jobs"]):
        raise RuntimeError("fresh observed job-count drift")

    out = {
        "schema": "ttf_genetic_phylogatr_phase3_shard_plan_v0.1",
        "status": "frozen_execution_plan",
        "geometry_fingerprint_sha256": auth["geometry_fingerprint_sha256"],
        "master_seed": int(auth["master_seed"]),
        "reference_jobs": reference_jobs,
        "observed_jobs": observed_jobs,
        "reference_job_count": len(reference_jobs),
        "observed_job_count": len(observed_jobs),
        "confirmatory_sequence_identity_opened": False,
        "confirmatory_pairwise_genetic_distances_opened": False,
        "qualification_claim_made": False,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(out, indent=2, sort_keys=True) + "\n")
    print(
        json.dumps(
            {
                "reference_job_count": len(reference_jobs),
                "observed_job_count": len(observed_jobs),
                "geometry_fingerprint_sha256": auth["geometry_fingerprint_sha256"],
            },
            sort_keys=True,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
