#!/usr/bin/env python3
from __future__ import annotations

import argparse,json
from pathlib import Path
import numpy as np

from ttf.private_null_inference import envelope_upper_pvalues
from ttf.lepidoptera_trait_gradient_v02 import wilson_interval

DESIGN_SHA="b1b9dab23676f6a1ea416c5889fc42676908a54d3abcf042dbaf8605cda2c3bf"
NUISANCE=(
    "private",
    "geometry_confounded_trap",
    "host_breadth_confounded_trap",
    "host_resource_gradient_trap",
)
CELLS=NUISANCE+("relational_latent_positive",)


def collect(directory:Path,namespace:str,cell:str,expected:int)->np.ndarray:
    rows=[]
    for path in sorted(directory.glob("*.json")):
        p=json.loads(path.read_text())
        if p.get("schema")!="ttf_lepidoptera_pair_repeatability_v01_qualification_shard":
            continue
        if p.get("namespace")!=namespace or p.get("cell")!=cell:
            continue
        if p.get("design_npz_sha256")!=DESIGN_SHA:
            raise RuntimeError("design hash drift")
        if int(p.get("retained_pairs",-1))!=11034:
            raise RuntimeError("retained pair count drift")
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

    references={cell:collect(args.input_dir,"reference",cell,499) for cell in NUISANCE}
    evaluation={}
    for cell in CELLS:
        values=collect(args.input_dir,"evaluation",cell,500)
        p,_=envelope_upper_pvalues(values,references)
        rejected=p<=.05
        k=int(np.count_nonzero(rejected))
        lo,hi=wilson_interval(k,len(values))
        evaluation[cell]={
            "worlds":500,
            "rejections":k,
            "rejection_rate":float(np.mean(rejected)),
            "wilson95_lower":lo,
            "wilson95_upper":hi,
            "mean_statistic":float(np.mean(values)),
            "median_statistic":float(np.median(values)),
            "median_envelope_p":float(np.median(p)),
            "minimum_envelope_p":float(np.min(p)),
        }

    gates={f"{cell}_type1_pass":evaluation[cell]["wilson95_upper"]<=.10 for cell in NUISANCE}
    gates["relational_positive_power_pass"]=evaluation["relational_latent_positive"]["wilson95_lower"]>=.80
    gates["all_five_pass"]=all(gates.values())

    out={
        "schema":"ttf_lepidoptera_pair_repeatability_v01_qualification_result",
        "status":"PASS" if gates["all_five_pass"] else "FAIL",
        "design_npz_sha256":DESIGN_SHA,
        "retained_pairs":11034,
        "reference":{
            cell:{
                "worlds":499,
                "mean_statistic":float(np.mean(values)),
                "median_statistic":float(np.median(values)),
            }
            for cell,values in references.items()
        },
        "evaluation":evaluation,
        "gates":gates,
        "empirical_split_half_repeatability_eligible":bool(gates["all_five_pass"]),
        "outcome_firewall":{
            "empirical_split_half_repeatability_computed":False,
        },
    }
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(out,indent=2,sort_keys=True)+"\n")
    print(json.dumps(gates,sort_keys=True))
    return 0


if __name__=="__main__":
    raise SystemExit(main())
