#!/usr/bin/env python3
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np

from ttf.ecological_information_ladder import ecological_information_ladder


def sha256_path(path: Path) -> str:
    h=hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda:f.read(1<<20),b""):
            h.update(chunk)
    return h.hexdigest()


def main()->int:
    ap=argparse.ArgumentParser()
    ap.add_argument("--historical-design",type=Path,required=True)
    ap.add_argument("--panel",choices=("development","confirmatory"),default="development")
    ap.add_argument("--output",type=Path,required=True)
    args=ap.parse_args()

    d=np.load(args.historical_design,allow_pickle=False)
    p=args.panel
    result=ecological_information_ladder(
        class_name=np.asarray(d["class_name"]),
        order_name=np.asarray(d["order_name"]),
        family_name=np.asarray(d["family_name"]),
        occurrence_n=np.asarray(d["retained_environment_occurrences"]),
        source_index=np.asarray(d[f"{p}_source_index"]),
        target_index=np.asarray(d[f"{p}_target_index"]),
        r_current=np.asarray(d[f"{p}_R_current"]),
        r_hist=np.asarray(d[f"{p}_R_hist"]),
    )
    payload={
        "status":"CHARACTERIZED_B1_TO_C1_RESPONSE_BLIND_INFORMATION_LADDER",
        "panel":p,
        "historical_design_sha256":sha256_path(args.historical_design),
        **result,
    }
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(payload,indent=2,sort_keys=True)+"\n")
    print(json.dumps({
        "status":payload["status"],
        "panel":p,
        "B1_total_unique_fraction":payload["layers"]["B1_current_climate"]["total_unique_variance_fraction"],
        "C1_total_unique_fraction":payload["layers"]["C1_history_given_current"]["total_unique_variance_fraction"],
        "hist_current_partial_correlation":payload["cross_layer"]["partial_correlation_hist_current_after_endpoint_lineage_locality"],
        "history_fraction_remaining_after_current":payload["cross_layer"]["history_fraction_remaining_after_adding_current_to_baseline"],
    },sort_keys=True))
    return 0


if __name__=="__main__":
    raise SystemExit(main())
