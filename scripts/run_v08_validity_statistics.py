#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

from run_geometry_calibration import columns, load_geometry
from ttf.batch import score_prepared_batch
from ttf.calibration import seed_for
from ttf.core import SpeciesSample, edge_turnover
from ttf.geometry import geometry_fingerprint, simulate_fixed_geometry_boundary_world
from ttf.geometry_control import length_orthogonalized_turnover
from ttf.nulls import edges_on_fixed_graphs, fixed_graphs
from ttf.private_strength import edge_midpoint_neighbor_indices, training_private_strength_from_indices
from ttf.transfer import prepare_transfer


def main() -> int:
    ap=argparse.ArgumentParser(description='Generate fresh v0.8 Gate-V private validity joint statistics.')
    ap.add_argument('--mode',choices=['observed','reference'],required=True)
    ap.add_argument('--panel',required=True)
    ap.add_argument('--input',type=Path,required=True)
    ap.add_argument('--source-ledger',type=Path,required=True)
    ap.add_argument('--coordinate-columns',type=columns,default=['x_km','y_km','z_km'])
    ap.add_argument('--min-records',type=int,default=80)
    ap.add_argument('--amplitude',type=float,required=True)
    ap.add_argument('--noise-sd',type=float,required=True)
    ap.add_argument('--configuration-label',default=None)
    ap.add_argument('--replicates',type=int,required=True)
    ap.add_argument('--k',type=int,default=4)
    ap.add_argument('--world-batch-size',type=int,default=25)
    ap.add_argument('--strength-neighbours',type=int,default=4)
    ap.add_argument('--transition-width',type=float,default=0.2)
    ap.add_argument('--output',type=Path,required=True)
    args=ap.parse_args()

    if args.mode=='reference' and not args.configuration_label:
        raise RuntimeError('reference mode requires configuration label')
    if args.mode=='observed' and args.configuration_label is not None:
        raise RuntimeError('observed mode must not have a configuration label')
    # Gate-V is deliberately type-I-only. A shared positive world is not accepted
    # by this runner; both modes always simulate shared_fraction=0 below.

    geometries,excluded=load_geometry(args.input,species_column='species',coordinate_columns=args.coordinate_columns,block_column=None,min_records=args.min_records)
    if excluded:
        raise RuntimeError(f'unexpected excluded taxa: {excluded}')
    gmap={g.species:g for g in geometries}
    ledger=json.loads(args.source_ledger.read_text())
    if ledger.get('candidate_performance_evaluated_at_freeze') is not False or ledger.get('synthetic_worlds_run_at_freeze') != 0:
        raise RuntimeError('validity geometry was not performance-unopened at freeze')
    exe=ledger['prospective_execution']
    train=tuple(exe['train_species'])
    evaluation=tuple(exe['eval_species'])
    if set(train)&set(evaluation) or set(train+evaluation)!=set(gmap):
        raise RuntimeError('frozen v0.8 split does not match geometry')
    if int(exe['k']) != args.k:
        raise RuntimeError('k drifted from frozen geometry contract')
    bandwidth=float(exe['bandwidth_km'])
    world_seed=int(exe['world_seed'])

    templates=[SpeciesSample(species=n,coordinates=gmap[n].coordinates,trait=np.arange(len(gmap[n].coordinates),dtype=float)) for n in train+evaluation]
    graphs=fixed_graphs(templates,k=args.k)
    template_edges=edges_on_fixed_graphs(templates,graphs)
    train_length={n:template_edges[n].length for n in train}
    strength_index={n:edge_midpoint_neighbor_indices(template_edges[n].midpoint,k=args.strength_neighbours) for n in train}
    prepared=prepare_transfer([template_edges[n] for n in train],[template_edges[n] for n in evaluation],bandwidth=bandwidth,prior_strength=0.25,prior_mean=0.0,segment_points=5)

    statistics=[]
    strengths=[]
    for start in range(0,args.replicates,args.world_batch_size):
        stop=min(start+args.world_batch_size,args.replicates)
        train_cols={n:[] for n in train}
        eval_cols={n:[] for n in evaluation}
        batch_strength=[]
        for replicate in range(start,stop):
            if args.mode=='observed':
                use_seed=seed_for(world_seed,'v08_validity_observed_private',float(args.amplitude),replicate)
            else:
                use_seed=seed_for(world_seed,'v08_validity_private_reference',args.configuration_label,replicate)
            world=simulate_fixed_geometry_boundary_world(geometries,shared_fraction=0.0,amplitude=float(args.amplitude),noise_sd=float(args.noise_sd),transition_width=float(args.transition_width),seed=use_seed)
            smap={s.species:s for s in world.samples}
            turnover={n:edge_turnover(smap[n],graphs[n]) for n in train+evaluation}
            train_response={n:length_orthogonalized_turnover(turnover[n],train_length[n]) for n in train}
            strength,_=training_private_strength_from_indices(train_response,strength_index,train)
            batch_strength.append(float(strength))
            for n in train: train_cols[n].append(train_response[n])
            for n in evaluation: eval_cols[n].append(turnover[n])
        scored=score_prepared_batch(prepared,{n:np.column_stack(train_cols[n]) for n in train},{n:np.column_stack(eval_cols[n]) for n in evaluation})
        statistics.extend(map(float,scored.statistics))
        strengths.extend(batch_strength)

    t=np.asarray(statistics,float); u=np.asarray(strengths,float)
    if t.shape!=(args.replicates,) or u.shape!=t.shape or not np.isfinite(t).all() or not np.isfinite(u).all():
        raise RuntimeError('incomplete/non-finite v0.8 validity statistics')

    payload={
        'schema':'ttf_v08_validity_statistics_v0.1','status':'fresh_external_validity_only',
        'mode':args.mode,'panel':args.panel,
        'geometry':{'input':str(args.input),'source_ledger':str(args.source_ledger),'fingerprint_sha256':geometry_fingerprint(geometries),'species_count':len(geometries),'record_count':int(sum(len(g.coordinates) for g in geometries))},
        'split':{'train_species':list(train),'eval_species':list(evaluation)},
        'config':{'shared_fraction':0.0,'amplitude':float(args.amplitude),'noise_sd':float(args.noise_sd),'configuration_label':args.configuration_label,'replicates':args.replicates,'k':args.k,'bandwidth':bandwidth,'world_seed':world_seed,'world_seed_namespace':'v08_validity_observed_private' if args.mode=='observed' else 'v08_validity_private_reference','strength_neighbours':args.strength_neighbours,'transition_width':args.transition_width},
        'statistics':statistics,'training_strength':strengths,
        'summary':{'statistic_mean':float(t.mean()),'statistic_sd':float(t.std(ddof=1)),'strength_mean':float(u.mean()),'strength_sd':float(u.std(ddof=1)),'strength_median':float(np.median(u)),'strength_mad':float(np.median(np.abs(u-np.median(u)))),'strength_q05':float(np.quantile(u,.05)),'strength_q95':float(np.quantile(u,.95)),'strength_statistic_correlation':float(np.corrcoef(u,t)[0,1])},
        'shared_positive_control_opened':False,'legacy_birds_butterflies_opened':False,'rgfca_reserve_opened':False,
        'claim_ready':False,
    }
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(payload,indent=2,sort_keys=True)+'\n')
    print(json.dumps({'panel':args.panel,'mode':args.mode,'label':args.configuration_label,'amplitude':args.amplitude,'summary':payload['summary']},sort_keys=True))
    return 0

if __name__=='__main__': raise SystemExit(main())
