#!/usr/bin/env python3
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np

from ttf.genetic_geometry import prepare_density_scaled_genetic_geometry
from ttf.lepidoptera_host_resource_qualification import (
    frozen_alignment_indices,
    geometry_kernel_from_features,
    residualize_host_resource_similarity,
    resource_breadth_kernel,
    resource_jaccard_kernel,
)

INPUT_SHA="597b1499dbdf8d4fb5a30f10a32686a60185305d12a8772019279bda08b82897"


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
    ap.add_argument("--census-design-npz",type=Path,required=True)
    ap.add_argument("--output-npz",type=Path,required=True)
    ap.add_argument("--output-json",type=Path,required=True)
    args=ap.parse_args()
    if sha256_path(args.census_design_npz)!=INPUT_SHA:
        raise RuntimeError("host-resource census design SHA drift")

    z=np.load(args.census_design_npz,allow_pickle=False)
    names=tuple(map(str,z["species_order"]))
    coordinates=np.asarray(z["coordinates"],float)
    offsets=np.asarray(z["coordinate_offsets"],np.int64)
    geos={}
    geometry_features=[]
    midpoint={}
    for i,name in enumerate(names):
        coords=coordinates[offsets[i]:offsets[i+1]]
        g=prepare_density_scaled_genetic_geometry(coords,neighbor_fraction=.15)
        geos[name]=g
        nodes=np.asarray(g.edge_nodes,np.int64)
        midpoint[name]=0.5*(coords[nodes[:,0]]+coords[nodes[:,1]])
        centroid=np.mean(coords,axis=0)
        extent=float(np.sqrt(np.sum(np.var(coords,axis=0))))
        geometry_features.append([
            *centroid,
            np.log1p(g.n_edges),
            np.log1p(g.n_localities),
            np.log1p(extent),
        ])
    geometry_features=np.asarray(geometry_features,float)
    geometry_kernel=geometry_kernel_from_features(geometry_features)

    host_presence=np.asarray(z["host_presence"],np.uint8)
    host_kernel=resource_jaccard_kernel(host_presence)
    native_hosts=np.asarray(z["hosts_with_primary_native_units"],float)
    footprint_units=np.sum(host_presence,axis=1,dtype=np.int64).astype(float)
    breadth_kernel=resource_breadth_kernel(native_hosts,footprint_units)

    target=np.asarray(z["target_index"],np.int64)
    source=np.asarray(z["source_index"],np.int64)
    raw=np.asarray(z["host_resource_jaccard"],float)
    if not np.allclose(raw,host_kernel[target,source],atol=1e-12,rtol=0):
        raise RuntimeError("host-resource Jaccard drift")

    nuisance=np.column_stack([
        np.asarray(z["coverage"],float),
        np.log1p(np.asarray(z["centroid_distance_km"],float)/500.0),
        np.log(np.asarray(z["edge_count_ratio"],float)),
        np.log(np.asarray(z["locality_count_ratio"],float)),
        geometry_kernel[target,source],
        np.log1p(native_hosts[source]),
        np.log1p(footprint_units[source]),
    ])
    residual=residualize_host_resource_similarity(raw,nuisance,target)

    align_target=[]
    align_source=[]
    align_offsets=[]
    align_lengths=[]
    cursor=0
    for ti,si in zip(target,source):
        tname=names[int(ti)]
        sname=names[int(si)]
        a,b=frozen_alignment_indices(
            midpoint[tname],
            midpoint[sname],
            target_name=tname,
            source_name=sname,
            radius=500.0,
            maximum_rows=128,
        )
        align_offsets.append(cursor)
        align_lengths.append(len(a))
        align_target.append(a)
        align_source.append(b)
        cursor+=len(a)
    flat_target=np.concatenate(align_target).astype(np.int64)
    flat_source=np.concatenate(align_source).astype(np.int64)

    per_target_sd=[]
    for label in np.unique(target):
        per_target_sd.append(float(np.std(residual[target==label],ddof=0)))

    args.output_npz.parent.mkdir(parents=True,exist_ok=True)
    np.savez_compressed(
        args.output_npz,
        species_order=np.asarray(names,dtype="U"),
        family=np.asarray(z["family"],dtype="U"),
        coordinates=coordinates,
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
        host_resource_kernel=host_kernel,
        resource_breadth_kernel=breadth_kernel,
        target_index=target,
        source_index=source,
        raw_host_resource_jaccard=raw,
        host_resource_residual=residual,
        nuisance_covariates=nuisance,
        alignment_offset=np.asarray(align_offsets,np.int64),
        alignment_length=np.asarray(align_lengths,np.int64),
        alignment_target=flat_target,
        alignment_source=flat_source,
    )

    out={
        "schema":"ttf_lepidoptera_host_resource_predictor_design_v0.3",
        "status":"FROZEN_RESPONSE_BLIND_BEFORE_CHARACTER_MASKS",
        "parent_census_design_sha256":INPUT_SHA,
        "output_design_npz_sha256":sha256_path(args.output_npz),
        "species":len(names),
        "pairs":len(target),
        "supported_eval_targets":len(np.unique(target)),
        "nuisance_covariates":[
            "target_to_source_coverage",
            "log1p_centroid_distance_over_500km",
            "log_edge_count_ratio",
            "log_locality_count_ratio",
            "pairwise_geometry_kernel_similarity",
            "log1p_source_native_host_species_count",
            "log1p_source_host_footprint_WGSRPD3_unit_count",
        ],
        "host_resource_predictor":{
            "raw_jaccard_quantiles":q(raw),
            "raw_jaccard_sd":float(np.std(raw)),
            "residual_sd":float(np.std(residual)),
            "residual_quantiles":q(residual),
            "within_target_residual_sd_quantiles":q(per_target_sd),
            "targets_with_residual_sd_gt_0_05":int(np.count_nonzero(np.asarray(per_target_sd)>.05)),
            "targets_with_residual_sd_gt_0_10":int(np.count_nonzero(np.asarray(per_target_sd)>.10)),
        },
        "alignment":{
            "radius_km":500.0,
            "maximum_rows_per_pair":128,
            "minimum_rows_per_pair":3,
            "rows_per_pair_quantiles":q(align_lengths),
            "total_alignment_rows":int(len(flat_target)),
            "deterministic_pair_specific_cap":True
        },
        "kernel_spread":{
            "geometry_pair_similarity_quantiles":q(geometry_kernel[target,source]),
            "resource_breadth_pair_similarity_quantiles":q(breadth_kernel[target,source]),
        },
        "outcome_firewall":{
            "character_masks_opened":False,
            "sequence_identity_opened":False,
            "pairwise_genetic_distances_opened":False,
            "transfer_statistic_computed":False
        }
    }
    args.output_json.parent.mkdir(parents=True,exist_ok=True)
    args.output_json.write_text(json.dumps(out,indent=2,sort_keys=True)+"\n")
    print(json.dumps({
        "design_sha256":out["output_design_npz_sha256"],
        "pairs":out["pairs"],
        "residual_sd":out["host_resource_predictor"]["residual_sd"],
        "targets_sd_gt_0_05":out["host_resource_predictor"]["targets_with_residual_sd_gt_0_05"],
        "alignment_rows":out["alignment"]["total_alignment_rows"]
    },sort_keys=True))
    return 0


if __name__=="__main__":
    raise SystemExit(main())
