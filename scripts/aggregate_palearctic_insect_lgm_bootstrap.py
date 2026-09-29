#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

from ttf.relational_qualification import relation_repeatability


def main() -> int:
    ap=argparse.ArgumentParser()
    ap.add_argument("--input-dir",type=Path,required=True)
    ap.add_argument("--output-matrix",type=Path,required=True)
    ap.add_argument("--output-summary",type=Path,required=True)
    args=ap.parse_args()

    rows=[]
    for replicate in range(100):
        matches=list(args.input_dir.rglob(f"relation_{replicate}.npy"))
        if len(matches)!=1:
            raise RuntimeError(
                f"expected exactly one bootstrap artifact for replicate {replicate}, found {len(matches)}"
            )
        rows.append(np.load(matches[0],allow_pickle=False))
    lengths={len(row) for row in rows}
    if len(lengths)!=1:
        raise RuntimeError("bootstrap relation dyad count drift")
    matrix=np.vstack(rows)
    summary=relation_repeatability(matrix)
    args.output_matrix.parent.mkdir(parents=True,exist_ok=True)
    np.save(args.output_matrix,matrix,allow_pickle=False)
    payload={
        "schema":"ttf_palearctic_insect_lgm_relation_repeatability_v0.1",
        "status":"COMPLETE_RESPONSE_BLIND_RELATION_REPEATABILITY",
        "replicates":100,
        "dyads":int(matrix.shape[1]),
        "repeatability":{
            "icc":float(summary.icc),
            "between_dyad_variance":float(summary.between_dyad_variance),
            "within_dyad_variance":float(summary.within_dyad_variance),
            "mean_within_dyad_sd":float(summary.mean_within_dyad_sd),
            "median_within_dyad_sd":float(summary.median_within_dyad_sd),
        },
        "genetic_response_used":False,
    }
    args.output_summary.write_text(json.dumps(payload,indent=2,sort_keys=True)+"\n")
    print(json.dumps(payload,sort_keys=True))
    return 0


if __name__=="__main__":
    raise SystemExit(main())
