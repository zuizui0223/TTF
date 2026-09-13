#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

from benchmark_genetic_decker_pipeline import load_geometries
from ttf.calibration import seed_for
from ttf.genetic_gate import prepare_genetic_ttf_design
from ttf.genetic_self_detectability import (
    prepare_genetic_self_detectability,
    score_genetic_self_world_batch,
)
from ttf.genetic_simulate import simulate_genetic_distance_world
from ttf.private_null_inference import upper_monte_carlo_pvalue


def main() -> int:
    ap=argparse.ArgumentParser()
    ap.add_argument('--geometry',type=Path,required=True)
    ap.add_argument('--manifest',type=Path,required=True)
    ap.add_argument('--pilot-rule',type=Path,required=True)
    ap.add_argument('--qualification',type=Path,required=True)
    ap.add_argument('--references',type=Path,required=True)
    ap.add_argument('--cell',choices=['null','private_A2'],required=True)
    ap.add_argument('--start',type=int,required=True)
    ap.add_argument('--stop',type=int,required=True)
    ap.add_argument('--output',type=Path,required=True)
    args=ap.parse_args()

    manifest=json.loads(args.manifest.read_text())
    pilot=json.loads(args.pilot_rule.read_text())
    qual=json.loads(args.qualification.read_text())
    refs=json.loads(args.references.read_text())
    if qual.get('schema')!='ttf_genetic_within_species_detectability_qualification_v0.1':
        raise RuntimeError('unexpected self qualification rule')
    if refs.get('schema')!='ttf_genetic_self_references_v0.1':
        raise RuntimeError('complete independent self reference required')
    if any(v is not False for v in qual['outcome_firewall'].values()):
        raise RuntimeError('empirical outcome firewall is open')
    syn=qual['synthetic_worlds']
    total=int(syn['null_evaluation_worlds'] if args.cell=='null' else syn['private_A2_evaluation_worlds'])
    start,stop=int(args.start),int(args.stop)
    if not 0 <= start < stop <= total:
        raise RuntimeError('evaluation range outside frozen universe')
    reference=np.asarray(refs['statistics'],dtype=float)
    if len(reference)!=int(syn['independent_null_reference_worlds']):
        raise RuntimeError('self reference length drift')

    geometries=load_geometries(args.geometry)
    split=manifest['split']
    design=prepare_genetic_ttf_design(
        geometries,
        train_species=split['train_species'],eval_species=split['eval_species'],
        bandwidth=float(pilot['core_v011_inheritance']['bandwidth_km']),
        prior_strength=float(pilot['core_v011_inheritance']['prior_strength']),
        segment_points=int(pilot['core_v011_inheritance']['segment_points']),
        min_training_edges=int(pilot['geometry_filter']['min_endpoint_disjoint_ibd_training_edges']),
        strength_neighbours=int(pilot['core_v011_inheritance']['strength_neighbours']),
    )
    self_design=prepare_genetic_self_detectability(
        design,bandwidth=float(qual['geometry']['bandwidth_km']),
        prior_strength=float(qual['geometry']['prior_strength']),
        prior_mean=float(qual['geometry']['prior_mean']),
        segment_points=int(qual['geometry']['segment_points']),
    )
    world_spec=syn['null_evaluation_world'] if args.cell=='null' else syn['private_A2_evaluation_world']
    seed_tag='genetic_self_null_eval' if args.cell=='null' else 'genetic_self_private_A2_eval'
    batch_width=int(qual['execution_only']['world_batch_width'])
    rows=[]
    for b0 in range(start,stop,batch_width):
        b1=min(b0+batch_width,stop)
        worlds=[
            simulate_genetic_distance_world(
                geometries,
                shared_fraction=float(world_spec['shared_fraction']),
                residual_amplitude=float(world_spec['residual_amplitude']),
                ibd_strength=float(syn['ibd_strength']),
                noise_sd=float(world_spec['noise_sd']),
                transition_width=float(syn['transition_width']),
                noise_dimensions=int(syn['noise_dimensions']),
                seed=seed_for(int(syn['master_seed']),seed_tag,rep),
            ) for rep in range(b0,b1)
        ]
        scored=score_genetic_self_world_batch(design,self_design,worlds)
        for col,rep in enumerate(range(b0,b1)):
            if int(scored.n_species[col])!=len(design.eval_species):
                raise RuntimeError('non-finite self species count')
            stat=float(scored.statistics[col])
            p=float(upper_monte_carlo_pvalue(stat,reference))
            rows.append({'replicate':int(rep),'statistic':stat,'p_value':p})
    if [r['replicate'] for r in rows]!=list(range(start,stop)):
        raise RuntimeError('self evaluation ordering drift')
    out={
        'schema':'ttf_genetic_self_evaluation_shard_v0.1',
        'status':'synthetic_self_detectability_evaluation_shard',
        'cell':args.cell,'start':start,'stop':stop,'n_worlds':stop-start,
        'alpha':float(qual['inference']['alpha']),'rows':rows,
        'empirical_genetic_outcomes_opened':False,'qualification_claim_made':False,
    }
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(out,indent=2,sort_keys=True)+'\n')
    print(json.dumps({k:v for k,v in out.items() if k!='rows'},sort_keys=True))
    return 0

if __name__=='__main__':
    raise SystemExit(main())
