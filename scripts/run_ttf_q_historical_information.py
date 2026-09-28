#!/usr/bin/env python3
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np

from ttf.historical_relation_information import historical_relation_information


def sha256_path(path: Path) -> str:
    h=hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda:handle.read(1<<20),b""):
            h.update(chunk)
    return h.hexdigest()


def main() -> int:
    ap=argparse.ArgumentParser()
    ap.add_argument("--historical-design",type=Path,required=True)
    ap.add_argument("--panel",choices=("development","confirmatory"),default="development")
    ap.add_argument("--output",type=Path,required=True)
    args=ap.parse_args()

    data=np.load(args.historical_design,allow_pickle=False)
    prefix=args.panel
    result=historical_relation_information(
        species_order=np.asarray(data["species_order"]),
        class_name=np.asarray(data["class_name"]),
        order_name=np.asarray(data["order_name"]),
        family_name=np.asarray(data["family_name"]),
        occurrence_n=np.asarray(data["retained_environment_occurrences"]),
        source_index=np.asarray(data[f"{prefix}_source_index"]),
        target_index=np.asarray(data[f"{prefix}_target_index"]),
        r_hist=np.asarray(data[f"{prefix}_R_hist"]),
        r_current=np.asarray(data[f"{prefix}_R_current"]),
    )
    payload={
        "schema":"ttf_q_historical_information_result_v0.1",
        "status":"CHARACTERIZED_HISTORY_SPECIFIC_ECOLOGICAL_INFORMATION",
        "panel":args.panel,
        "historical_design_sha256":sha256_path(args.historical_design),
        "program_boundary":{
            "relational_v03_reopened":False,
            "genetic_response_used":False,
            "geographic_opportunity_used":False,
            "role":"fresh response-blind ecological characterization",
        },
        **result,
        "interpretation":{
            "high_current_low_history":"present-day niche convergence reached through dissimilar late-Quaternary displacement histories",
            "high_current_high_history":"present-day similarity accompanied by similar late-Quaternary displacement histories",
            "nonredundancy":"history-specific pairwise information after current climate, lineage, locality imbalance, and source/target identity",
            "not_yet_full_ttf_q":"geographic-opportunity adjustment belongs to the later full qualification stage",
        },
    }
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(payload,indent=2,sort_keys=True)+"\n")
    print(json.dumps({
        "status":payload["status"],
        "panel":args.panel,
        "spearman_hist_vs_current":payload["ecological_geometry"]["spearman_hist_vs_current"],
        "high_current_low_history_fraction":payload["ecological_geometry"]["high_current_low_history_fraction"],
        "total_unique_variance_fraction":payload["nonredundancy"]["total_unique_variance_fraction"],
        "signal_effective_sources":payload["dyadic_signal_support"]["signal_effective_sources"],
        "signal_effective_targets":payload["dyadic_signal_support"]["signal_effective_targets"],
    },sort_keys=True))
    return 0


if __name__=="__main__":
    raise SystemExit(main())
