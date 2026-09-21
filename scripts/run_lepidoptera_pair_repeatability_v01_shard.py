#!/usr/bin/env python3
from __future__ import annotations

import argparse,hashlib,json
from pathlib import Path
import numpy as np

from ttf.genetic_geometry import prepare_density_scaled_genetic_geometry
from ttf.lepidoptera_pair_repeatability import (
    prepare_pair_repeatability,
    score_pair_repeatability,
)
from ttf.lepidoptera_pair_repeatability_simulate import (
    ALLOWED_CELLS,
    EVALUATION_TAG,
    REFERENCE_TAG,
    frozen_seed,
    prepare_pair_repeatability_simulator,
    simulate_pair_repeatability_world,
)

DESIGN_SHA="b1b9dab23676f6a1ea416c5889fc42676908a54d3abcf042dbaf8605cda2c3bf"
NUISANCE_CELLS=(
    "private",
    "geometry_confounded_trap",
    "host_breadth_confounded_trap",
    "host_resource_gradient_trap",
)


def sha256_path(path:Path)->str:
    h=hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda:f.read(1<<20),b""):
            h.update(chunk)
    return h.hexdigest()


def load_design(path:Path):
    if sha256_path(path)!=DESIGN_SHA:
        raise RuntimeError("pair-repeatability design SHA drift")
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
    edge_counts=np.asarray([geos[n].n_edges for n in names],np.int64)
    repeatability=prepare_pair_repeatability(
        z,edge_counts,minimum_rows_per_half=8
    )
    if len(repeatability.pair_index)!=11034:
        raise RuntimeError("retained pair count drift")
    simulator=prepare_pair_repeatability_simulator(
        geos,names,
        geometry_kernel=np.asarray(z["geometry_kernel"],float),
        resource_breadth_kernel=np.asarray(z["resource_breadth_kernel"],float),
        host_resource_kernel=np.asarray(z["host_resource_positive_kernel"],float),
        relational_dimensions=8,
    )
    return z,names,repeatability,simulator


def main()->int:
    ap=argparse.ArgumentParser()
    ap.add_argument("--design-npz",type=Path,required=True)
    ap.add_argument("--namespace",choices=("reference","evaluation"),required=True)
    ap.add_argument("--cell",choices=ALLOWED_CELLS,required=True)
    ap.add_argument("--start",type=int,required=True)
    ap.add_argument("--count",type=int,required=True)
    ap.add_argument("--output",type=Path,required=True)
    args=ap.parse_args()

    if args.start<0 or args.count<1 or args.count>500:
        raise RuntimeError("formal shards require start>=0 and 1<=count<=500")
    if args.namespace=="reference" and args.cell not in NUISANCE_CELLS:
        raise RuntimeError("relational positive has no reference distribution")
    if args.namespace=="reference" and args.start+args.count>499:
        raise RuntimeError("reference shard exceeds 499 worlds")
    if args.namespace=="evaluation" and args.start+args.count>500:
        raise RuntimeError("evaluation shard exceeds 500 worlds")

    _,names,repeatability,simulator=load_design(args.design_npz)
    tag=REFERENCE_TAG if args.namespace=="reference" else EVALUATION_TAG
    stats=[]; fractions=[]; covariances=[]
    for replicate in range(args.start,args.start+args.count):
        world=simulate_pair_repeatability_world(
            simulator,
            cell=args.cell,
            seed=frozen_seed(20260921,tag,args.cell,replicate),
            latent_fields=6,
            private_amplitude=.35,
            noise_sd=.10,
            transition_width=.20,
        )
        score=score_pair_repeatability(
            world.edge_response,names,repeatability
        )
        stats.append(float(score.statistic))
        fractions.append(float(score.pair_variance_fraction))
        covariances.append(float(score.covariance))

    payload={
        "schema":"ttf_lepidoptera_pair_repeatability_v01_qualification_shard",
        "namespace":args.namespace,
        "cell":args.cell,
        "start":args.start,
        "count":args.count,
        "design_npz_sha256":DESIGN_SHA,
        "retained_pairs":int(len(repeatability.pair_index)),
        "statistics":stats,
        "pair_variance_fractions":fractions,
        "covariances":covariances,
        "outcome_firewall":{
            "empirical_split_half_repeatability_computed":False,
        },
    }
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(payload,indent=2)+"\n")
    print(json.dumps({
        "namespace":args.namespace,
        "cell":args.cell,
        "start":args.start,
        "count":args.count,
        "retained_pairs":payload["retained_pairs"],
    },sort_keys=True))
    return 0


if __name__=="__main__":
    raise SystemExit(main())
