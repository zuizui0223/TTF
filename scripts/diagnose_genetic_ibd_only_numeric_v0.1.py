#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

from benchmark_genetic_decker_pipeline import load_geometries
from ttf.core import average_ranks
from ttf.genetic_gate import prepare_genetic_ttf_design, score_genetic_world
from ttf.genetic_ibd import crossfit_ibd_residuals
from ttf.genetic_simulate import simulate_genetic_distance_world


def main() -> int:
    ap=argparse.ArgumentParser()
    ap.add_argument('--geometry',type=Path,required=True)
    ap.add_argument('--manifest',type=Path,required=True)
    ap.add_argument('--rule',type=Path,required=True)
    ap.add_argument('--authorization',type=Path,required=True)
    ap.add_argument('--output',type=Path,required=True)
    args=ap.parse_args()

    manifest=json.loads(args.manifest.read_text())
    rule=json.loads(args.rule.read_text())
    auth=json.loads(args.authorization.read_text())
    if auth.get('schema')!='ttf_genetic_ibd_only_numeric_diagnostic_authorization_v0.1':
        raise RuntimeError('authorization drift')
    world_rule=auth['authorized_world']
    geometries=load_geometries(args.geometry)
    world=simulate_genetic_distance_world(
        geometries,
        shared_fraction=float(world_rule['shared_fraction']),
        residual_amplitude=float(world_rule['residual_amplitude']),
        ibd_strength=float(world_rule['ibd_strength']),
        noise_sd=float(world_rule['noise_sd']),
        transition_width=float(rule['genetic_world']['transition_width']),
        noise_dimensions=int(rule['genetic_world']['noise_dimensions']),
        seed=int(world_rule['seed']),
    )

    per_species={}
    all_rank_identical=True
    max_abs=[]
    ranges=[]
    unique_counts=[]
    for name in sorted(geometries):
        g=geometries[name]
        nodes=g.edge_nodes
        geo=np.linalg.norm(g.coordinates[nodes[:,0]]-g.coordinates[nodes[:,1]],axis=1)
        genetic=np.asarray(world.genetic_distance[name],dtype=float)
        rank_identical=bool(np.array_equal(average_ranks(genetic),average_ranks(geo)))
        result=crossfit_ibd_residuals(
            genetic,geo,nodes,
            min_training_edges=int(rule['geometry_filter']['min_endpoint_disjoint_ibd_training_edges']),
        )
        ma=float(np.max(np.abs(result.residual)))
        rr=float(np.ptp(result.residual))
        uc=int(len(np.unique(result.residual_turnover)))
        all_rank_identical &= rank_identical
        max_abs.append(ma); ranges.append(rr); unique_counts.append(uc)
        per_species[name]={
            'global_rank_order_identical':rank_identical,
            'max_abs_raw_residual':ma,
            'raw_residual_range':rr,
            'unique_post_rank_turnover':uc,
            'n_edges':int(len(nodes)),
        }

    design=prepare_genetic_ttf_design(
        geometries,
        train_species=manifest['split']['train_species'],
        eval_species=manifest['split']['eval_species'],
        bandwidth=float(rule['core_v011_inheritance']['bandwidth_km']),
        prior_strength=float(rule['core_v011_inheritance']['prior_strength']),
        segment_points=int(rule['core_v011_inheritance']['segment_points']),
        min_training_edges=int(rule['geometry_filter']['min_endpoint_disjoint_ibd_training_edges']),
        strength_neighbours=int(rule['core_v011_inheritance']['strength_neighbours']),
    )
    score=score_genetic_world(design,world)
    payload={
        'schema':'ttf_genetic_ibd_only_numeric_diagnostic_v0.1',
        'status':'ibd_only_diagnostic_not_qualification',
        'species_count':len(geometries),
        'all_species_global_rank_order_identical':all_rank_identical,
        'species_with_rank_order_mismatch':[
            n for n,r in per_species.items() if not r['global_rank_order_identical']
        ],
        'raw_residual':{
            'max_of_species_max_abs':float(max(max_abs)),
            'median_species_max_abs':float(np.median(max_abs)),
            'max_species_range':float(max(ranges)),
            'median_species_range':float(np.median(ranges)),
        },
        'post_rank_turnover':{
            'median_unique_values_per_species':float(np.median(unique_counts)),
            'max_unique_values_per_species':int(max(unique_counts)),
            'species_with_more_than_one_unique_value':int(np.count_nonzero(np.asarray(unique_counts)>1)),
        },
        'current_pipeline':{
            'heldout_statistic':float(score.statistic),
            'training_strength':float(score.training_strength),
        },
        'per_species':per_species,
        'private_or_shared_performance_opened':False,
        'empirical_genetic_outcomes_opened':False,
        'threshold_grid_search_performed':False,
    }
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(payload,indent=2,sort_keys=True)+'\n')
    print(json.dumps({k:v for k,v in payload.items() if k!='per_species'},sort_keys=True))
    return 0

if __name__=='__main__':
    raise SystemExit(main())
