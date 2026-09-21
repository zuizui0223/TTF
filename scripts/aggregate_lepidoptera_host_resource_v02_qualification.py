#!/usr/bin/env python3
from __future__ import annotations

import argparse,json
from pathlib import Path
import numpy as np

from ttf.lepidoptera_host_resource_qualification_v02 import nuisance_envelope_pvalues
from ttf.lepidoptera_trait_gradient_v02 import wilson_interval

DESIGN_SHA="b1b9dab23676f6a1ea416c5889fc42676908a54d3abcf042dbaf8605cda2c3bf"
CELLS=("private","geometry_confounded_trap","host_breadth_confounded_trap","host_resource_gradient_positive")


def collect(directory:Path,namespace:str,cell:str,expected:int):
    rows=[]
    for path in sorted(directory.glob("*.json")):
        p=json.loads(path.read_text())
        if p.get("schema")!="ttf_lepidoptera_host_resource_v02_qualification_shard":
            continue
        if p.get("namespace")!=namespace or p.get("cell")!=cell:
            continue
        if p.get("design_npz_sha256")!=DESIGN_SHA:
            raise RuntimeError("v0.2 design hash drift")
        start=int(p["start"])
        vals=list(map(float,p["statistics"]))
        if len(vals)!=int(p["count"]):
            raise RuntimeError("shard count drift")
        rows.extend((start+i,v) for i,v in enumerate(vals))
    if len(rows)!=expected:
        raise RuntimeError(f"{namespace}/{cell}: {len(rows)} != {expected}")
    rows.sort()
    if [i for i,_ in rows]!=list(range(expected)):
        raise RuntimeError(f"{namespace}/{cell}: replicate coverage drift")
    return np.asarray([v for _,v in rows],float)


def main()->int:
    ap=argparse.ArgumentParser()
    ap.add_argument("--input-dir",type=Path,required=True)
    ap.add_argument("--output",type=Path,required=True)
    args=ap.parse_args()

    rp=collect(args.input_dir,"reference","private",499)
    rg=collect(args.input_dir,"reference","geometry_confounded_trap",499)
    rb=collect(args.input_dir,"reference","host_breadth_confounded_trap",499)
    evaluation={}
    for cell in CELLS:
        v=collect(args.input_dir,"evaluation",cell,500)
        p=nuisance_envelope_pvalues(
            v,private_reference=rp,geometry_reference=rg,breadth_reference=rb
        )
        rejected=p<=.05
        k=int(np.count_nonzero(rejected))
        lo,hi=wilson_interval(k,len(v))
        evaluation[cell]={
            "worlds":500,
            "rejections":k,
            "rejection_rate":float(np.mean(rejected)),
            "wilson95_lower":lo,
            "wilson95_upper":hi,
            "mean_statistic":float(np.mean(v)),
            "median_statistic":float(np.median(v)),
            "median_envelope_p":float(np.median(p)),
            "minimum_envelope_p":float(np.min(p)),
        }

    gates={
        "private_type1_pass":evaluation["private"]["wilson95_upper"]<=.10,
        "geometry_trap_type1_pass":evaluation["geometry_confounded_trap"]["wilson95_upper"]<=.10,
        "host_breadth_trap_type1_pass":evaluation["host_breadth_confounded_trap"]["wilson95_upper"]<=.10,
        "host_resource_positive_power_pass":evaluation["host_resource_gradient_positive"]["wilson95_lower"]>=.80,
    }
    gates["all_four_pass"]=all(gates.values())
    out={
        "schema":"ttf_lepidoptera_host_resource_v02_qualification_result",
        "status":"PASS" if gates["all_four_pass"] else "FAIL",
        "design_npz_sha256":DESIGN_SHA,
        "reference":{
            "private_worlds":499,
            "geometry_confounded_trap_worlds":499,
            "host_breadth_confounded_trap_worlds":499,
            "private_mean_statistic":float(np.mean(rp)),
            "geometry_mean_statistic":float(np.mean(rg)),
            "host_breadth_mean_statistic":float(np.mean(rb)),
        },
        "evaluation":evaluation,
        "gates":gates,
        "empirical_identity_opening_eligible":bool(gates["all_four_pass"]),
        "outcome_firewall":{
            "sequence_identity_opened":False,
            "pairwise_genetic_distances_opened":False,
            "transfer_statistic_computed":False,
        }
    }
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(out,indent=2,sort_keys=True)+"\n")
    print(json.dumps(gates,sort_keys=True))
    return 0


if __name__=="__main__":
    raise SystemExit(main())
