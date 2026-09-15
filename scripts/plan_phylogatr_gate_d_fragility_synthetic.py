#!/usr/bin/env python3
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

from ttf.phylogatr_fragility_execution import build_fragility_execution_plan
from ttf.phylogatr_fragility_formal_anchor import validate_formal_phase3_qualification


EXECUTION_RULE_GIT_BLOB_SHA = "a06eac22f47d90a29328a837c401a6c8b938c6ab"


def _git_blob_sha1(path: Path) -> str:
    data = path.read_bytes()
    header = f"blob {len(data)}\0".encode("ascii")
    return hashlib.sha1(header + data).hexdigest()


def _apply_identical_geometry_policy(plan: dict) -> dict:
    """Prevent Monte Carlo noise from creating fake fragility at geometry plateaus."""
    seen: dict[str, float] = {
        str(plan["full_geometry_fingerprint_sha256"]): 1.0,
    }
    for level in plan["levels"]:
        retention = float(level["retention_fraction"])
        fingerprint = str(level["geometry_fingerprint_sha256"])
        inherited_from = seen.get(fingerprint)
        if inherited_from is None:
            seen[fingerprint] = retention
            continue
        level["status"] = (
            "IDENTICAL_TO_FULL_GEOMETRY_NO_RERUN"
            if inherited_from == 1.0
            else "IDENTICAL_TO_HIGHER_RETENTION_GEOMETRY_NO_RERUN"
        )
        level["inherits_metrics_from_retention_fraction"] = inherited_from
        level["reference_jobs"] = []
        level["observed_jobs"] = []
        level["reference_job_count"] = 0
        level["observed_job_count"] = 0
    return plan


def main() -> int:
    ap = argparse.ArgumentParser(
        description="Plan diagnostic-only synthetic Gate-D fragility execution for thinned geometries."
    )
    ap.add_argument("--geometry-plan-receipt", type=Path, required=True)
    ap.add_argument("--phase2-manifest", type=Path, required=True)
    ap.add_argument("--phase3-rule", type=Path, required=True)
    ap.add_argument("--formal-qualification", type=Path, required=True)
    ap.add_argument("--execution-rule", type=Path, required=True)
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()

    if _git_blob_sha1(args.execution_rule) != EXECUTION_RULE_GIT_BLOB_SHA:
        raise RuntimeError(
            "frozen fragility execution rule Git blob SHA drift before diagnostic planning"
        )

    phase2 = json.loads(args.phase2_manifest.read_text())
    expected_fingerprint = str(phase2["geometry_fingerprint_sha256"])
    validate_formal_phase3_qualification(
        args.phase3_rule,
        args.formal_qualification,
        expected_geometry_fingerprint_sha256=expected_fingerprint,
    )

    plan = build_fragility_execution_plan(
        args.geometry_plan_receipt,
        args.phase2_manifest,
        args.phase3_rule,
        args.formal_qualification,
        args.execution_rule,
    )
    plan = _apply_identical_geometry_policy(plan)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(plan, indent=2, sort_keys=True) + "\n")
    print(
        json.dumps(
            {
                "status": plan["status"],
                "formal_anchor_status": plan["formal_anchor"]["formal_status"],
                "levels": [
                    {
                        "retention_fraction": level["retention_fraction"],
                        "status": level["status"],
                        "inherits_metrics_from_retention_fraction": level.get(
                            "inherits_metrics_from_retention_fraction"
                        ),
                        "reference_job_count": level["reference_job_count"],
                        "observed_job_count": level["observed_job_count"],
                    }
                    for level in plan["levels"]
                ],
                "formal_gate_d_decision_made_by_this_plan": False,
                "phase4_identity_opening_authorized_by_this_plan": False,
            },
            sort_keys=True,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
