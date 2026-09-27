#!/usr/bin/env python3
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np

from ttf.historical_geographic_information import (
    nested_historical_geographic_information,
)


def sha256_path(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--opportunity-design", type=Path, required=True)
    ap.add_argument("--c1-result", type=Path, required=True)
    ap.add_argument(
        "--rule",
        type=Path,
        default=Path("docs/supporting/ttf_q_c2_geographic_opportunity_v0.1.json"),
    )
    ap.add_argument(
        "--panel",
        choices=("development", "confirmatory"),
        required=True,
    )
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()

    rule = json.loads(args.rule.read_text())
    c1 = json.loads(args.c1_result.read_text())
    if rule.get("schema") != "ttf_q_c2_geographic_opportunity_v0.1":
        raise RuntimeError("unexpected TTF-Q C2 rule")
    if c1.get("schema") != "ttf_q_historical_information_result_v0.1":
        raise RuntimeError("unexpected C1 result")
    if c1.get("panel") != args.panel:
        raise RuntimeError("C1 panel does not match requested C2 panel")
    if c1["program_boundary"]["genetic_response_used"]:
        raise RuntimeError("C1 response firewall is open")

    data = np.load(args.opportunity_design, allow_pickle=False)
    p = args.panel
    result = nested_historical_geographic_information(
        source_index=np.asarray(data[f"{p}_source_index"]),
        target_index=np.asarray(data[f"{p}_target_index"]),
        r_hist=np.asarray(data[f"{p}_R_hist"]),
        r_current=np.asarray(data[f"{p}_R_current"]),
        geographic_coverage=np.asarray(data[f"{p}_coverage"]),
        same_class=np.asarray(data[f"{p}_same_class"]),
        same_order=np.asarray(data[f"{p}_same_order"]),
        same_family=np.asarray(data[f"{p}_same_family"]),
        occurrence_count_ratio=np.asarray(data[f"{p}_locality_count_ratio"]),
    )

    recomputed = result["C1_recomputed"]["total_unique_variance_fraction"]
    observed = c1["nonredundancy"]["total_unique_variance_fraction"]
    if not np.isclose(recomputed, observed, rtol=1e-10, atol=1e-12):
        raise RuntimeError(
            f"C1 nesting drift: recomputed={recomputed} frozen={observed}"
        )

    payload = {
        "schema": "ttf_q_c2_geographic_information_result_v0.1",
        "status": "CHARACTERIZED_GEOGRAPHY_SPECIFIC_INFORMATION_LOSS",
        "panel": args.panel,
        "opportunity_design_sha256": sha256_path(args.opportunity_design),
        "c1_result_sha256": sha256_path(args.c1_result),
        "rule_sha256": sha256_path(args.rule),
        "program_boundary": {
            "relational_v03_reopened": False,
            "genetic_response_used": False,
            "role": "fresh response-blind ecological characterization",
        },
        **result,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n"
    )
    print(json.dumps({
        "status": payload["status"],
        "panel": args.panel,
        "C1_total_unique": observed,
        "C2_total_unique": payload["C2_geography_adjusted"][
            "total_unique_variance_fraction"
        ],
        "geography_retained_fraction": payload["nested_information"][
            "geography_retained_fraction"
        ],
        "geography_attributable_fraction": payload["nested_information"][
            "geography_attributable_fraction_of_C1_unique_information"
        ],
    }, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
