#!/usr/bin/env python3
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np

from ttf.genetic_geometry import prepare_density_scaled_genetic_geometry
from ttf.lepidoptera_trait_gradient_simulate import (
    prepare_trait_gradient_simulator,
    simulate_prepared_trait_gradient_world,
)
from ttf.lepidoptera_trait_gradient_v02 import (
    EVALUATION_TAG,
    REFERENCE_TAG,
    frozen_v02_seed,
    target_equal_mean_correlation,
    target_fe_geometry_residual,
)


def sha256_path(path: Path) -> str:
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def load_design(path: Path):
    z=np.load(path,allow_pickle=False)
    names=tuple(map(str,z["species_order"]))
    coords=np.asarray(z["coordinates"],float)
    offsets=np.asarray(z["coordinate_offsets"],np.int64)
    pairs=np.asarray(z["pairs"],float)
    target=pairs[:,0].astype(np.int64)
    source=pairs[:,1].astype(np.int64)
    geometry=np.asarray(pairs[:,3:7],float).copy()
    geometry[:,1]=np.log1p(geometry[:,1]/500.0)
    geometry[:,2]=np.log(geometry[:,2])
    geometry[:,3]=np.log(geometry[:,3])
    geometry=np.column_stack([
        geometry,
        np.asarray(z["geometry_kernel"],float)[target,source],
    ])
    residual=target_fe_geometry_residual(pairs[:,2],geometry,target)
    geos={
        name:prepare_density_scaled_genetic_geometry(
            coords[offsets[i]:offsets[i+1]],neighbor_fraction=0.15
        )
        for i,name in enumerate(names)
    }
    edge_counts=np.asarray([geos[name].n_edges for name in names],np.int64)
    edge_offsets=np.r_[0,np.cumsum(edge_counts)]
    lengths=pairs[:,8].astype(np.int64)
    pair_id=np.repeat(np.arange(len(pairs),dtype=np.int64),lengths)
    expected_offsets=np.r_[0,np.cumsum(lengths)[:-1]]
    if not np.array_equal(pairs[:,7].astype(np.int64),expected_offsets):
        raise RuntimeError("alignment offset drift")
    align_target=np.asarray(z["alignment_target"],np.int64)
    align_source=np.asarray(z["alignment_source"],np.int64)
    if len(align_target)!=len(pair_id) or len(align_source)!=len(pair_id):
        raise RuntimeError("alignment length drift")
    global_target=edge_offsets[target][pair_id]+align_target
    global_source=edge_offsets[source][pair_id]+align_source
    return z,names,geos,pairs,target,residual,pair_id,global_target,global_source


def pair_transfer_congruence(world,names,pair_id,global_target,global_source,n_pairs):
    flat=np.concatenate([np.asarray(world.edge_response[name],float) for name in names])
    x=flat[global_target]; y=flat[global_source]
    n=np.bincount(pair_id,minlength=n_pairs).astype(float)
    sx=np.bincount(pair_id,weights=x,minlength=n_pairs)
    sy=np.bincount(pair_id,weights=y,minlength=n_pairs)
    sxx=np.bincount(pair_id,weights=x*x,minlength=n_pairs)
    syy=np.bincount(pair_id,weights=y*y,minlength=n_pairs)
    sxy=np.bincount(pair_id,weights=x*y,minlength=n_pairs)
    covariance=sxy-sx*sy/n
    denominator=np.sqrt(np.maximum((sxx-sx*sx/n)*(syy-sy*sy/n),0.0))
    correlation=np.zeros(n_pairs,dtype=float)
    valid=denominator>np.sqrt(np.finfo(float).eps)
    correlation[valid]=covariance[valid]/denominator[valid]
    return correlation


def main() -> int:
    ap=argparse.ArgumentParser()
    ap.add_argument("--design-npz",type=Path,required=True)
    ap.add_argument("--namespace",choices=["reference","evaluation"],required=True)
    ap.add_argument("--cell",choices=["private","geometry_confounded_trap","trait_gradient_positive"],required=True)
    ap.add_argument("--start",type=int,required=True)
    ap.add_argument("--count",type=int,required=True)
    ap.add_argument("--output",type=Path,required=True)
    args=ap.parse_args()
    if args.start<0 or args.count<1 or args.count>50:
        raise RuntimeError("v0.2 shards require start>=0 and 1<=count<=50")
    if args.namespace=="reference" and args.cell=="trait_gradient_positive":
        raise RuntimeError("positive cell has no reference distribution")
    if args.namespace=="reference" and args.start+args.count>499:
        raise RuntimeError("reference shard exceeds 499 worlds")
    if args.namespace=="evaluation" and args.start+args.count>500:
        raise RuntimeError("evaluation shard exceeds 500 worlds")

    z,names,geos,pairs,target,residual,pair_id,global_target,global_source=load_design(args.design_npz)
    tag=REFERENCE_TAG if args.namespace=="reference" else EVALUATION_TAG
    fraction={"private":0.0,"geometry_confounded_trap":0.85,"trait_gradient_positive":0.90}[args.cell]
    statistics=[]
    target_counts=[]
    for replicate in range(args.start,args.start+args.count):
        world=simulate_trait_gradient_world(
            geos,
            names,
            np.asarray(z["trait_kernel"],float),
            np.asarray(z["geometry_kernel"],float),
            cell=args.cell,
            seed=frozen_v02_seed(20260920,tag,args.cell,replicate),
            shared_fraction=fraction,
            private_amplitude=0.35,
            noise_sd=0.10,
            transition_width=0.20,
            latent_fields=6,
        )
        pair_score=pair_transfer_congruence(
            world,names,pair_id,global_target,global_source,len(pairs)
        )
        statistic,per_target=target_equal_mean_correlation(pair_score,residual,target)
        statistics.append(statistic)
        target_counts.append(len(per_target))

    payload={
        "schema":"ttf_lepidoptera_trait_gradient_v02_shard",
        "namespace":args.namespace,
        "cell":args.cell,
        "start":args.start,
        "count":args.count,
        "design_npz_sha256":sha256_path(args.design_npz),
        "statistics":list(map(float,statistics)),
        "target_count_min":int(min(target_counts)),
        "target_count_max":int(max(target_counts)),
        "outcome_firewall":{
            "new_lepidoptera_empirical_sequence_identity_opened":False,
            "new_lepidoptera_empirical_pairwise_genetic_distances_opened":False,
            "new_lepidoptera_empirical_transfer_statistic_computed":False,
        },
    }
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(payload,indent=2)+"\n")
    print(json.dumps({
        "namespace":args.namespace,"cell":args.cell,"start":args.start,
        "count":args.count,"mean_statistic":float(np.mean(statistics)),
        "target_count_min":payload["target_count_min"],
    },sort_keys=True))
    return 0


if __name__=="__main__":
    raise SystemExit(main())
