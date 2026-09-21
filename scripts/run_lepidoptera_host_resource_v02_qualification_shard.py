#!/usr/bin/env python3
from __future__ import annotations

import argparse,hashlib,json
from pathlib import Path
import numpy as np

from ttf.genetic_geometry import prepare_density_scaled_genetic_geometry
from ttf.lepidoptera_host_resource_qualification_v02 import (
    ALLOWED_CELLS,
    EVALUATION_TAG,
    REFERENCE_TAG,
    frozen_v02_seed,
)
from ttf.lepidoptera_host_resource_simulate import (
    prepare_host_resource_simulator,
    simulate_prepared_host_resource_world,
)
from ttf.lepidoptera_trait_gradient_v02 import target_equal_mean_correlation

DESIGN_SHA="b1b9dab23676f6a1ea416c5889fc42676908a54d3abcf042dbaf8605cda2c3bf"


def sha256_path(path:Path)->str:
    h=hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda:f.read(1<<20),b""):
            h.update(chunk)
    return h.hexdigest()


def load_design(path:Path):
    if sha256_path(path)!=DESIGN_SHA:
        raise RuntimeError("v0.2 host-resource design SHA drift")
    z=np.load(path,allow_pickle=False)
    names=tuple(map(str,z["species_order"]))
    coords=np.asarray(z["coordinates"],float)
    offsets=np.asarray(z["coordinate_offsets"],np.int64)
    geos={
        name:prepare_density_scaled_genetic_geometry(
            coords[offsets[i]:offsets[i+1]],neighbor_fraction=.15
        )
        for i,name in enumerate(names)
    }
    target=np.asarray(z["target_index"],np.int64)
    source=np.asarray(z["source_index"],np.int64)
    residual=np.asarray(z["host_resource_residual"],float)
    lengths=np.asarray(z["alignment_length"],np.int64)
    pair_offsets=np.asarray(z["alignment_offset"],np.int64)
    expected=np.r_[0,np.cumsum(lengths)[:-1]]
    if not np.array_equal(pair_offsets,expected):
        raise RuntimeError("alignment offset drift")
    pair_id=np.repeat(np.arange(len(target),dtype=np.int64),lengths)
    align_target=np.asarray(z["alignment_target"],np.int64)
    align_source=np.asarray(z["alignment_source"],np.int64)
    edge_counts=np.asarray([geos[n].n_edges for n in names],np.int64)
    edge_offsets=np.r_[0,np.cumsum(edge_counts)]
    global_target=edge_offsets[target][pair_id]+align_target
    global_source=edge_offsets[source][pair_id]+align_source
    simulator=prepare_host_resource_simulator(
        geos,names,
        geometry_kernel=np.asarray(z["geometry_kernel"],float),
        resource_breadth_kernel=np.asarray(z["resource_breadth_kernel"],float),
        host_resource_positive_kernel=np.asarray(z["host_resource_positive_kernel"],float),
    )
    return z,names,simulator,target,residual,pair_id,global_target,global_source


def pair_transfer_congruence(world,names,pair_id,global_target,global_source,n_pairs):
    flat=np.concatenate([np.asarray(world.edge_response[n],float) for n in names])
    x=flat[global_target]; y=flat[global_source]
    n=np.bincount(pair_id,minlength=n_pairs).astype(float)
    sx=np.bincount(pair_id,weights=x,minlength=n_pairs)
    sy=np.bincount(pair_id,weights=y,minlength=n_pairs)
    sxx=np.bincount(pair_id,weights=x*x,minlength=n_pairs)
    syy=np.bincount(pair_id,weights=y*y,minlength=n_pairs)
    sxy=np.bincount(pair_id,weights=x*y,minlength=n_pairs)
    cov=sxy-sx*sy/n
    den=np.sqrt(np.maximum((sxx-sx*sx/n)*(syy-sy*sy/n),0.0))
    out=np.zeros(n_pairs,float)
    ok=den>np.sqrt(np.finfo(float).eps)
    out[ok]=cov[ok]/den[ok]
    return out


def main()->int:
    ap=argparse.ArgumentParser()
    ap.add_argument("--design-npz",type=Path,required=True)
    ap.add_argument("--namespace",choices=("reference","evaluation"),required=True)
    ap.add_argument("--cell",choices=ALLOWED_CELLS,required=True)
    ap.add_argument("--start",type=int,required=True)
    ap.add_argument("--count",type=int,required=True)
    ap.add_argument("--output",type=Path,required=True)
    args=ap.parse_args()
    if args.start<0 or args.count<1 or args.count>100:
        raise RuntimeError("v0.2 formal shards require start>=0 and 1<=count<=100")
    if args.namespace=="reference" and args.cell=="host_resource_gradient_positive":
        raise RuntimeError("positive cell has no nuisance reference")
    if args.namespace=="reference" and args.start+args.count>499:
        raise RuntimeError("reference shard exceeds 499 worlds")
    if args.namespace=="evaluation" and args.start+args.count>500:
        raise RuntimeError("evaluation shard exceeds 500 worlds")

    z,names,simulator,target,residual,pair_id,gt,gs=load_design(args.design_npz)
    tag=REFERENCE_TAG if args.namespace=="reference" else EVALUATION_TAG
    statistics=[]
    counts=[]
    for replicate in range(args.start,args.start+args.count):
        world=simulate_prepared_host_resource_world(
            simulator,
            cell=args.cell,
            seed=frozen_v02_seed(20260921,tag,args.cell,replicate),
            latent_fields=6,
            private_amplitude=.35,
            noise_sd=.10,
            transition_width=.20,
        )
        pair_score=pair_transfer_congruence(world,names,pair_id,gt,gs,len(target))
        statistic,per_target=target_equal_mean_correlation(pair_score,residual,target)
        statistics.append(float(statistic))
        counts.append(len(per_target))

    payload={
        "schema":"ttf_lepidoptera_host_resource_v02_qualification_shard",
        "namespace":args.namespace,
        "cell":args.cell,
        "start":args.start,
        "count":args.count,
        "design_npz_sha256":DESIGN_SHA,
        "statistics":statistics,
        "target_count_min":int(min(counts)),
        "target_count_max":int(max(counts)),
        "outcome_firewall":{
            "sequence_identity_opened":False,
            "pairwise_genetic_distances_opened":False,
            "transfer_statistic_computed":False,
        }
    }
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(payload,indent=2)+"\n")
    print(json.dumps({
        "namespace":args.namespace,
        "cell":args.cell,
        "start":args.start,
        "count":args.count,
        "target_count_min":payload["target_count_min"],
    },sort_keys=True))
    return 0


if __name__=="__main__":
    raise SystemExit(main())
