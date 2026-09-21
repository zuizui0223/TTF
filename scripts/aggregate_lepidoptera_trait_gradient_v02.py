#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

from ttf.lepidoptera_trait_gradient_v02 import envelope_pvalues,wilson_interval


CELLS=("private","geometry_confounded_trap","trait_gradient_positive")


def collect(directory:Path,namespace:str,cell:str,expected:int):
    rows=[]
    design_sha=None
    for path in sorted(directory.glob("*.json")):
        payload=json.loads(path.read_text())
        if payload.get("schema")!="ttf_lepidoptera_trait_gradient_v02_shard":
            continue
        if payload.get("namespace")!=namespace or payload.get("cell")!=cell:
            continue
        if design_sha is None:
            design_sha=payload["design_npz_sha256"]
        elif design_sha!=payload["design_npz_sha256"]:
            raise RuntimeError("design hash drift across shards")
        start=int(payload["start"]); values=list(map(float,payload["statistics"]))
        if len(values)!=int(payload["count"]):
            raise RuntimeError("shard count drift")
        rows.extend((start+i,value) for i,value in enumerate(values))
    if len(rows)!=expected:
        raise RuntimeError(f"{namespace}/{cell}: {len(rows)} != {expected}")
    rows.sort()
    if [i for i,_ in rows]!=list(range(expected)):
        raise RuntimeError(f"{namespace}/{cell}: replicate coverage drift")
    return np.asarray([value for _,value in rows],float),design_sha


def main()->int:
    ap=argparse.ArgumentParser()
    ap.add_argument("--input-dir",type=Path,required=True)
    ap.add_argument("--output",type=Path,required=True)
    args=ap.parse_args()

    ref_private,sha=collect(args.input_dir,"reference","private",499)
    ref_geometry,sha2=collect(args.input_dir,"reference","geometry_confounded_trap",499)
    if sha!=sha2:
        raise RuntimeError("reference design hash drift")

    evaluations={}
    for cell in CELLS:
        values,cell_sha=collect(args.input_dir,"evaluation",cell,500)
        if cell_sha!=sha:
            raise RuntimeError("evaluation design hash drift")
        p=envelope_pvalues(values,ref_private,ref_geometry)
        rejected=p<=0.05
        k=int(np.count_nonzero(rejected))
        lo,hi=wilson_interval(k,len(values))
        evaluations[cell]={
            "worlds":len(values),
            "rejections":k,
            "rejection_rate":float(np.mean(rejected)),
            "wilson95_lower":lo,
            "wilson95_upper":hi,
            "mean_statistic":float(np.mean(values)),
            "median_statistic":float(np.median(values)),
            "median_envelope_p":float(np.median(p)),
            "minimum_envelope_p":float(np.min(p)),
        }

    private_pass=evaluations["private"]["wilson95_upper"]<=0.10
    trap_pass=evaluations["geometry_confounded_trap"]["wilson95_upper"]<=0.10
    positive_pass=evaluations["trait_gradient_positive"]["wilson95_lower"]>=0.80
    passed=private_pass and trap_pass and positive_pass
    payload={
        "schema":"ttf_lepidoptera_trait_gradient_v02_qualification_result",
        "status":"PASS" if passed else "FAIL",
        "design_npz_sha256":sha,
        "reference":{
            "private_worlds":len(ref_private),
            "geometry_confounded_trap_worlds":len(ref_geometry),
            "private_mean_statistic":float(np.mean(ref_private)),
            "geometry_mean_statistic":float(np.mean(ref_geometry)),
        },
        "evaluation":evaluations,
        "gates":{
            "private_type1_pass":private_pass,
            "geometry_trap_type1_pass":trap_pass,
            "trait_gradient_positive_power_pass":positive_pass,
            "all_three_pass":passed,
        },
        "empirical_opening_authorized":bool(passed),
        "outcome_firewall":{
            "new_lepidoptera_empirical_sequence_identity_opened":False,
            "new_lepidoptera_empirical_pairwise_genetic_distances_opened":False,
            "new_lepidoptera_empirical_transfer_statistic_computed":False,
        },
    }
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(payload,indent=2,sort_keys=True)+"\n")
    print(json.dumps(payload["gates"],sort_keys=True))
    return 0


if __name__=="__main__":
    raise SystemExit(main())
