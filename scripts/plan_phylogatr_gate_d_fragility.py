#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
from pathlib import Path

from ttf.phylogatr_fragility import build_fragility_plan, write_fragility_plan


def main() -> int:
    ap = argparse.ArgumentParser(
        description=(
            "Create response-blind, diagnostic-only nested locality-thinning geometries "
            "for fresh phylogatR Gate-D fragility analysis."
        )
    )
    ap.add_argument("--geometry", type=Path, required=True)
    ap.add_argument("--phase2-manifest", type=Path, required=True)
    ap.add_argument("--phase3-rule", type=Path, required=True)
    ap.add_argument("--fragility-rule", type=Path, required=True)
    ap.add_argument("--output-dir", type=Path, required=True)
    args = ap.parse_args()

    plan = build_fragility_plan(
        args.geometry,
        args.phase2_manifest,
        args.phase3_rule,
        args.fragility_rule,
    )
    receipt_path = write_fragility_plan(plan, args.output_dir)
    receipt = json.loads(receipt_path.read_text())
    print(
        json.dumps(
            {
                "status": receipt["status"],
                "levels": [
                    {
                        "retention_fraction": level["retention_fraction"],
                        "minimum_localities": level["metrics"]["minimum_localities"],
                        "minimum_endpoint_disjoint_ibd_training_edges": level["metrics"][
                            "minimum_endpoint_disjoint_ibd_training_edges"
                        ],
                    }
                    for level in receipt["levels"]
                ],
                "formal_gate_d_decision_made_by_this_receipt": False,
                "phase4_identity_opening_authorized_by_this_receipt": False,
            },
            sort_keys=True,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
