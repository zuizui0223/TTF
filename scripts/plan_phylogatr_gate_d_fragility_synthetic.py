#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
from pathlib import Path

from ttf.phylogatr_fragility_execution import build_fragility_execution_plan


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

    plan = build_fragility_execution_plan(
        args.geometry_plan_receipt,
        args.phase2_manifest,
        args.phase3_rule,
        args.formal_qualification,
        args.execution_rule,
    )
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
