#!/usr/bin/env python3
from __future__ import annotations

import argparse,hashlib,json
from pathlib import Path
import numpy as np

from ttf.genetic_geometry import prepare_density_scaled_genetic_geometry
from ttf.lepidoptera_host_resource_qualification import (
    frozen_alignment_indices,
    geometry_kernel_from_features,
    residualize_host_resource_similarity,
    resource_breadth_kernel,
    resource_cosine_kernel,
    resource_jaccard_kernel,
)

INPUT_SHA="cef1fb5d34476692a0168eaf527ccea8b676b68409b88cdcf40db202d62f6e6c"


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
    ap.add_argument("--survivor-design-npz",type=Path,required=True)
    ap.add_argument("--output-npz",type=Path,required=True)
    ap.add_argument("--output-json",type=Path,required=True)
    args=ap.parse_args()
    if sha256_path(args.survivor_design_npz)!=INPUT_SHA:
        raise RuntimeError("survivor design SHA drift")

    z=np.load(args.survivor_design_npz,allow_pickle=False)
    names=tuple(map(str,z["species_order"]))
    coords=np.asarray(z["coordinates"],float)
    offsets=np.asarray(z["coordinate_offsets"],np.int64)
    host_presence=np.asarray(z["host_presence"],np.uint8)
    native_hosts=np.asarray(z["hosts_with_primary_native_units"],float)
    target=np.asarray(z["support_target_index"],np.int64)
    source=np.asarray(z["support_source_index"],np.int64)
    coverage=np.asarray(z["support_coverage"],float)

    geos={}
    midpoint={}
    geometry_features=[]
    for i,name in enumerate(names):
        x=coords[offsets[i]:offsets[i+1]]
        g=prepare_density_scaled_genetic_geometry(x,neighbor_fraction=.15)
        geos[name]=g
        nodes=np.asarray(g.edge_nodes,np.int64)
        midpoint[name]=0.5*(x[nodes[:,0]]+x[nodes[:,1]])
        centroid=np.mean(x,axis=0)
        extent=float(np.sqrt(np.sum(np.var(x,axis=0))))
        geometry_features.append([
            *centroid,
            np.log1p(g.n_edges),
            np.log1p(g.n_localities),
            np.log1p(extent),
        ])
    geometry_features=np.asarray(geometry_features,float)
    geometry_kernel=geometry_kernel_from_features(geometry_features)
    jaccard_kernel=resource_jaccard_kernel(host_presence)
    cosine_kernel=resource_cosine_kernel(host_presence)
    footprint_units=np.sum(host_presence,axis=1,dtype=np.int64).astype(float)
    breadth_kernel=resource_breadth_kernel(native_hosts,footprint_units)

    raw=jaccard_kernel[target,source]
    centroid_distance=[]
    edge_ratio=[]
    locality_ratio=[]
    for ti,si in zip(target,source):
        tg=geos[names[int(ti)]]
        sg=geos[names[int(si)]]
        centroid_distance.append(float(np.linalg.norm(
            np.mean(tg.coordinates,axis=0)-np.mean(sg.coordinates,axis=0)
        )))
        edge_ratio.append(float(sg.n_edges/tg.n_edges))
        locality_ratio.append(float(sg.n_localities/tg.n_localities))
    centroid_distance=np.asarray(centroid_distance,float)
    edge_ratio=np.asarray(edge_ratio,float)
    locality_ratio=np.asarray(locality_ratio,float)

    nuisance=np.column_stack([
        coverage,
        np.log1p(centroid_distance/500.0),
        np.log(edge_ratio),
        np.log(locality_ratio),
        geometry_kernel[target,source],
        np.log1p(native_hosts[source]),
        np.log1p(footprint_units[source]),
    ])
    residual=residualize_host_resource_similarity(raw,nuisance,target)

    align_target=[]
    align_source=[]
    align_offset=[]
    align_length=[]
    cursor=0
    for ti,si in zip(target,source):
        tn=names[int(ti)]
        sn=names[int(si)]
        a,b=frozen_alignment_indices(
            midpoint[tn],midpoint[sn],
            target_name=tn,source_name=sn,
            radius=500.0,maximum_rows=128,
        )
        align_offset.append(cursor)
        align_length.append(len(a))
        align_target.append(a)
        align_source.append(b)
        cursor+=len(a)
    flat_target=np.concatenate(align_target).astype(np.int64)
    flat_source=np.concatenate(align_source).astype(np.int64)

    per_target_sd=[
        float(np.std(residual[target==label],ddof=0))
        for label in np.unique(target)
    ]

    args.output_npz.parent.mkdir(parents=True,exist_ok=True)
    np.savez_compressed(
        args.output_npz,
        species_order=np.asarray(names,dtype="U"),
        family=np.asarray(z["family"],dtype="U"),
        coordinates=coords,
        coordinate_offsets=offsets,
        train_indices=np.asarray(z["train_indices"],np.int64),
        eval_indices=np.asarray(z["eval_indices"],np.int64),
        eligible_eval_indices=np.asarray(z["eligible_eval_indices"],np.int64),
        host_unit_names=np.asarray(z["host_unit_names"],dtype="U"),
        host_presence=host_presence,
        native_host_species_count=native_hosts.astype(np.int64),
        footprint_unit_count=footprint_units.astype(np.int64),
        geometry_features=geometry_features,
        geometry_kernel=geometry_kernel,
        host_resource_jaccard_kernel=jaccard_kernel,
        host_resource_positive_kernel=cosine_kernel,
        resource_breadth_kernel=breadth_kernel,
        target_index=target,
        source_index=source,
        coverage=coverage,
        centroid_distance_km=centroid_distance,
        edge_count_ratio=edge_ratio,
        locality_count_ratio=locality_ratio,
        raw_host_resource_jaccard=raw,
        host_resource_residual=residual,
        nuisance_covariates=nuisance,
        alignment_offset=np.asarray(align_offset,np.int64),
        alignment_length=np.asarray(align_length,np.int64),
        alignment_target=flat_target,
        alignment_source=flat_source,
    )

    out={
        "schema":"ttf_lepidoptera_host_resource_survivor_predictor_design_v0.1",
        "status":"FROZEN_RESPONSE_BLIND_EXACT_SURVIVOR_DESIGN",
        "parent_survivor_design_sha256":INPUT_SHA,
        "output_design_npz_sha256":sha256_path(args.output_npz),
        "species":len(names),
        "train_species":len(np.asarray(z["train_indices"])),
        "eval_species":len(np.asarray(z["eval_indices"])),
        "supported_eval_targets":len(np.unique(target)),
        "pairs":len(target),
        "primary_predictor":"raw native host-resource WGSRPD3 Jaccard, geometry/breadth-residualized within target",
        "positive_world_kernel":"PSD cosine similarity of binary native host-resource WGSRPD3 presence",
        "host_resource_predictor":{
            "raw_jaccard_quantiles":q(raw),
            "raw_jaccard_sd":float(np.std(raw)),
            "residual_quantiles":q(residual),
            "residual_sd":float(np.std(residual)),
            "within_target_residual_sd_quantiles":q(per_target_sd),
            "targets_with_residual_sd_gt_0_05":int(np.count_nonzero(np.asarray(per_target_sd)>.05)),
            "targets_with_residual_sd_gt_0_10":int(np.count_nonzero(np.asarray(per_target_sd)>.10)),
        },
        "alignment":{
            "rows_per_pair_quantiles":q(align_length),
            "total_alignment_rows":int(len(flat_target)),
            "maximum_rows_per_pair":128,
            "minimum_rows_per_pair":3,
        },
        "kernel_spread":{
            "geometry_pair_similarity_quantiles":q(geometry_kernel[target,source]),
            "breadth_pair_similarity_quantiles":q(breadth_kernel[target,source]),
            "positive_cosine_pair_similarity_quantiles":q(cosine_kernel[target,source]),
        },
        "outcome_firewall":{
            "sequence_identity_opened":False,
            "pairwise_genetic_distances_opened":False,
            "transfer_statistic_computed":False,
        }
    }
    args.output_json.parent.mkdir(parents=True,exist_ok=True)
    args.output_json.write_text(json.dumps(out,indent=2,sort_keys=True)+"\n")
    print(json.dumps({
        "design_sha":out["output_design_npz_sha256"],
        "species":out["species"],
        "pairs":out["pairs"],
        "supported_eval_targets":out["supported_eval_targets"],
        "residual_sd":out["host_resource_predictor"]["residual_sd"],
        "alignment_rows":out["alignment"]["total_alignment_rows"]
    },sort_keys=True))
    return 0


if __name__=="__main__":
    raise SystemExit(main())
