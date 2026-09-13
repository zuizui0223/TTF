#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
from pathlib import Path

from benchmark_genetic_decker_pipeline import load_geometries
from ttf.calibration import seed_for
from ttf.genetic_gate import prepare_genetic_ttf_design
from ttf.genetic_self_detectability import (
    prepare_genetic_self_detectability,
    score_genetic_self_world_batch,
)
from ttf.genetic_simulate import simulate_genetic_distance_world


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument('--geometry', type=Path, required=True)
    ap.add_argument('--manifest', type=Path, required=True)
    ap.add_argument('--pilot-rule', type=Path, required=True)
    ap.add_argument('--qualification', type=Path, required=True)
    ap.add_argument('--start', type=int, required=True)
    ap.add_argument('--stop', type=int, required=True)
    ap.add_argument('--output', type=Path, required=True)
    args = ap.parse_args()

    manifest=json.loads(args.manifest.read_text())
    pilot=json.loads(args.pilot_rule.read_text())
    qual=json.loads(args.qualification.read_text())
    if qual.get('schema')!='ttf_genetic_within_species_detectability_qualification_v0.1':
        raise RuntimeError('unexpected self-detectability qualification rule')
    if any(v is not False for v in qual['outcome_firewall'].values()):
        raise RuntimeError('empirical outcome firewall is open')
    total=int(qual['synthetic_worlds']['independent_null_reference_worlds'])
    start,stop=int(args.start),int(args.stop)
    if not 0 <= start < stop <= total:
        raise RuntimeError('reference range outside frozen universe')

    geometries=load_geometries(args.geometry)
    split=manifest['split']
    design=prepare_genetic_ttf_design(
        geometries,
        train_species=split['train_species'],
        eval_species=split['eval_species'],
        bandwidth=float(pilot['core_v011_inheritance']['bandwidth_km']),
        prior_strength=float(pilot['core_v011_inheritance']['prior_strength']),
        segment_points=int(pilot['core_v011_inheritance']['segment_points']),
        min_training_edges=int(pilot['geometry_filter']['min_endpoint_disjoint_ibd_training_edges']),
        strength_neighbours=int(pilot['core_v011_inheritance']['strength_neighbours']),
    )
    self_design=prepare_genetic_self_detectability(
        design,
        bandwidth=float(qual['geometry']['bandwidth_km']),
        prior_strength=float(qual['geometry']['prior_strength']),
        prior_mean=float(qual['geometry']['prior_mean']),
        segment_points=int(qual['geometry']['segment_points']),
    )
    batch_width=int(qual['execution_only']['world_batch_width'])
    syn=qual['synthetic_worlds']
    rows=[]
    for b0 in range(start,stop,batch_width):
        b1=min(b0+batch_width,stop)
        worlds=[
            simulate_genetic_distance_world(
                geometries,
                shared_fraction=float(syn['reference_world']['shared_fraction']),
                residual_amplitude=float(syn['reference_world']['residual_amplitude']),
                ibd_strength=float(syn['ibd_strength']),
                noise_sd=float(syn['reference_world']['noise_sd']),
                transition_width=float(syn['transition_width']),
                noise_dimensions=int(syn['noise_dimensions']),
                seed=seed_for(int(syn['master_seed']),'genetic_self_reference',rep),
            )
            for rep in range(b0,b1)
        ]
        scored=score_genetic_self_world_batch(design,self_design,worlds)
        for col,rep in enumerate(range(b0,b1)):
            if int(scored.n_species[col]) != len(design.eval_species):
                raise RuntimeError('non-finite species count in self reference')
            rows.append({'replicate':int(rep),'statistic':float(scored.statistics[col])})
    if [r['replicate'] for r in rows] != list(range(start,stop)):
        raise RuntimeError('reference ordering drift')
    out={
        'schema':'ttf_genetic_self_reference_shard_v0.1',
        'status':'independent_null_reference_shard',
        'start':start,'stop':stop,'n_worlds':stop-start,'rows':rows,
        'empirical_genetic_outcomes_opened':False,
        'qualification_claim_made':False,
    }
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(out,indent=2,sort_keys=True)+'\n')
    print(json.dumps({k:v for k,v in out.items() if k!='rows'},sort_keys=True))
    return 0

if __name__=='__main__':
    raise SystemExit(main())
