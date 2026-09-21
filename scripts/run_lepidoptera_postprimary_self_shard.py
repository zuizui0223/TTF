#!/usr/bin/env python3
from __future__ import annotations

import argparse, hashlib, json
from pathlib import Path
from types import SimpleNamespace

import numpy as np

from ttf.calibration import seed_for
from ttf.core import SpeciesEdges
from ttf.genetic_geometry import prepare_density_scaled_genetic_geometry
from ttf.genetic_self_detectability import prepare_genetic_self_detectability
from ttf.genetic_simulate import simulate_genetic_distance_world
from ttf.phylogatr_compact_ibd import prepare_compact_crossfit_ibd_design
from ttf.phylogatr_compact_self_detectability import score_phylogatr_compact_self_world_batch
from ttf.private_null_inference import upper_monte_carlo_pvalue

RULE_SCHEMA="ttf_lepidoptera_postprimary_self_detectability_rule_v0.1"


def sha256_path(path: Path) -> str:
    h=hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda:f.read(1<<20),b""):
            h.update(chunk)
    return h.hexdigest()


def template_edges(name: str, geometry):
    nodes=np.asarray(geometry.edge_nodes,np.int64)
    start=geometry.coordinates[nodes[:,0]]
    end=geometry.coordinates[nodes[:,1]]
    return SpeciesEdges(
        species=name,
        nodes=nodes,
        start=start,
        end=end,
        midpoint=.5*(start+end),
        length=np.linalg.norm(end-start,axis=1),
        turnover=np.zeros(len(nodes),float),
    )


def prepare_self_from_npz(path: Path):
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
    eval_idx=np.asarray(z["eval_indices"],np.int64)
    eval_names=tuple(names[int(i)] for i in eval_idx)
    templates={name:template_edges(name,geos[name]) for name in eval_names}
    ibd={
        name:prepare_compact_crossfit_ibd_design(
            templates[name].length,
            templates[name].nodes,
            min_training_edges=5,
        )
        for name in eval_names
    }
    design=SimpleNamespace(
        template_edges=templates,
        ibd_designs=ibd,
        eval_species=eval_names,
        min_training_edges=5,
    )
    self_design=prepare_genetic_self_detectability(
        design,
        species=eval_names,
        bandwidth=500.0,
        prior_strength=.25,
        prior_mean=.5,
        segment_points=5,
    )
    return z,names,geos,eval_names,design,self_design


def main()->int:
    ap=argparse.ArgumentParser()
    ap.add_argument("--design-npz",type=Path,required=True)
    ap.add_argument("--rule",type=Path,required=True)
    ap.add_argument("--panel",choices=("butterfly","host"),required=True)
    ap.add_argument("--cell",choices=("reference","null","private_A2"),required=True)
    ap.add_argument("--start",type=int,required=True)
    ap.add_argument("--count",type=int,required=True)
    ap.add_argument("--reference-npy",type=Path)
    ap.add_argument("--output",type=Path,required=True)
    args=ap.parse_args()

    rule=json.loads(args.rule.read_text())
    if rule.get("schema")!=RULE_SCHEMA:
        raise RuntimeError("self-detectability rule schema drift")
    panel_key="butterfly_trait_panel" if args.panel=="butterfly" else "host_resource_panel"
    panel_rule=rule["panels"][panel_key]
    if sha256_path(args.design_npz)!=panel_rule["survivor_design_sha256"]:
        raise RuntimeError("panel design SHA drift")

    n_total=1000 if args.cell=="reference" else 500
    if not 0<=args.start<args.start+args.count<=n_total:
        raise RuntimeError("self-detectability shard range drift")

    z,names,geos,eval_names,design,self_design=prepare_self_from_npz(args.design_npz)
    syn=rule["exact_geometry_qualification_per_panel"]
    spec=(
        syn["reference_world"] if args.cell=="reference"
        else syn["null_evaluation_world"] if args.cell=="null"
        else syn["private_A2_evaluation_world"]
    )
    tag=rule["seed_namespaces"][
        f"{args.panel}_{'reference' if args.cell=='reference' else 'null' if args.cell=='null' else 'private_A2'}"
    ]
    reference=None
    if args.cell!="reference":
        if args.reference_npy is None:
            raise RuntimeError("evaluation cell requires reference")
        reference=np.load(args.reference_npy,allow_pickle=False)
        if reference.shape!=(1000,):
            raise RuntimeError("self reference length drift")

    batch_width=32
    stats=[]
    pvals=[]
    for batch_start in range(args.start,args.start+args.count,batch_width):
        batch_stop=min(batch_start+batch_width,args.start+args.count)
        worlds=[
            simulate_genetic_distance_world(
                geos,
                shared_fraction=float(spec["shared_fraction"]),
                residual_amplitude=float(spec["residual_amplitude"]),
                ibd_strength=float(syn["ibd_strength"]),
                noise_sd=float(spec["noise_sd"]),
                transition_width=float(syn["transition_width"]),
                noise_dimensions=int(syn["noise_dimensions"]),
                seed=seed_for(20260921,tag,replicate),
            )
            for replicate in range(batch_start,batch_stop)
        ]
        scored=score_phylogatr_compact_self_world_batch(design,self_design,worlds)
        if np.any(scored.n_species!=len(eval_names)):
            raise RuntimeError("non-finite self species count")
        for value in scored.statistics:
            stats.append(float(value))
            if reference is not None:
                pvals.append(float(upper_monte_carlo_pvalue(float(value),reference)))

    out={
        "schema":"ttf_lepidoptera_postprimary_self_shard_v0.1",
        "panel":args.panel,
        "cell":args.cell,
        "start":args.start,
        "count":args.count,
        "design_npz_sha256":sha256_path(args.design_npz),
        "evaluation_species":len(eval_names),
        "statistics":stats,
        "p_values":pvals if reference is not None else None,
        "cross_species_primary_rewritten":False,
    }
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(out,indent=2)+"\n")
    print(json.dumps({
        "panel":args.panel,"cell":args.cell,"start":args.start,
        "count":args.count,"evaluation_species":len(eval_names)
    },sort_keys=True))
    return 0


if __name__=="__main__":
    raise SystemExit(main())
