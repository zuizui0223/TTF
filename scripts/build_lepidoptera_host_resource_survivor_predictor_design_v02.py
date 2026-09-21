#!/usr/bin/env python3
from __future__ import annotations

import argparse,hashlib,json
from pathlib import Path
import numpy as np

from ttf.lepidoptera_trait_gradient_v02 import target_fe_geometry_residual

INPUT_SHA="60cadfab5b8cf57462ee0e03be92ba44de52c628ad66e2bc8e209f1a40d786a2"


def sha256_path(path:Path)->str:
    h=hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda:f.read(1<<20),b""):
            h.update(chunk)
    return h.hexdigest()


def q(values):
    x=np.asarray(values,float)
    z=np.quantile(x,[0,.25,.5,.75,1])
    return dict(zip(("min","q25","median","q75","max"),map(float,z)))


def main()->int:
    ap=argparse.ArgumentParser()
    ap.add_argument("--v01-design-npz",type=Path,required=True)
    ap.add_argument("--output-npz",type=Path,required=True)
    ap.add_argument("--output-json",type=Path,required=True)
    args=ap.parse_args()
    if sha256_path(args.v01_design_npz)!=INPUT_SHA:
        raise RuntimeError("v0.1 design SHA drift")

    z=np.load(args.v01_design_npz,allow_pickle=False)
    target=np.asarray(z["target_index"],np.int64)
    source=np.asarray(z["source_index"],np.int64)
    raw=np.asarray(z["raw_host_resource_jaccard"],float)
    old_nuisance=np.asarray(z["nuisance_covariates"],float)
    breadth_pair=np.asarray(z["resource_breadth_kernel"],float)[target,source]
    nuisance=np.column_stack([old_nuisance,breadth_pair])
    residual=target_fe_geometry_residual(raw,nuisance,target)

    data={key:z[key] for key in z.files}
    data["nuisance_covariates"]=nuisance
    data["host_resource_residual"]=residual
    args.output_npz.parent.mkdir(parents=True,exist_ok=True)
    np.savez_compressed(args.output_npz,**data)

    per_target_sd=np.asarray([
        np.std(residual[target==label],ddof=0)
        for label in np.unique(target)
    ],float)
    breadth_center=breadth_pair-float(np.mean(breadth_pair))
    residual_center=residual-float(np.mean(residual))
    corr=float(
        np.dot(breadth_center,residual_center)
        /np.sqrt(np.dot(breadth_center,breadth_center)*np.dot(residual_center,residual_center))
    )
    out={
        "schema":"ttf_lepidoptera_host_resource_survivor_predictor_design_v0.2",
        "status":"FROZEN_RESPONSE_BLIND_V02_SUCCESSOR_DESIGN",
        "parent_v01_design_sha256":INPUT_SHA,
        "output_design_npz_sha256":sha256_path(args.output_npz),
        "single_change":"append pairwise resource-breadth similarity to the frozen nuisance covariates before target-fixed-effect residualization",
        "pairs":int(len(target)),
        "supported_eval_targets":int(len(np.unique(target))),
        "residual":{
            "sd":float(np.std(residual)),
            "quantiles":q(residual),
            "correlation_with_pairwise_resource_breadth_similarity":corr,
            "within_target_sd_quantiles":q(per_target_sd),
            "targets_with_sd_gt_0_05":int(np.count_nonzero(per_target_sd>.05)),
            "targets_with_sd_gt_0_10":int(np.count_nonzero(per_target_sd>.10))
        },
        "outcome_firewall":{
            "sequence_identity_opened":False,
            "pairwise_genetic_distances_opened":False,
            "transfer_statistic_computed":False
        }
    }
    args.output_json.parent.mkdir(parents=True,exist_ok=True)
    args.output_json.write_text(json.dumps(out,indent=2,sort_keys=True)+"\n")
    print(json.dumps({
        "design_sha":out["output_design_npz_sha256"],
        "residual_sd":out["residual"]["sd"],
        "breadth_corr":out["residual"]["correlation_with_pairwise_resource_breadth_similarity"]
    },sort_keys=True))
    return 0


if __name__=="__main__":
    raise SystemExit(main())
