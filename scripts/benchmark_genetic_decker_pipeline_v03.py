#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import time
from pathlib import Path

import numpy as np

from benchmark_genetic_decker_pipeline import load_geometries, sha256_path
from ttf.genetic_gate import prepare_genetic_ttf_design, score_genetic_world
from ttf.genetic_simulate import simulate_genetic_distance_world
from ttf.geometry import SpeciesGeometry, geometry_fingerprint


def main() -> int:
    ap=argparse.ArgumentParser()
    ap.add_argument('--geometry',type=Path,required=True)
    ap.add_argument('--manifest',type=Path,required=True)
    ap.add_argument('--source-ledger',type=Path,required=True)
    ap.add_argument('--rule',type=Path,required=True)
    ap.add_argument('--successor-rule',type=Path,required=True)
    ap.add_argument('--authorization',type=Path,required=True)
    ap.add_argument('--output',type=Path,required=True)
    args=ap.parse_args()

    manifest=json.loads(args.manifest.read_text())
    ledger=json.loads(args.source_ledger.read_text())
    rule=json.loads(args.rule.read_text())
    successor=json.loads(args.successor_rule.read_text())
    auth=json.loads(args.authorization.read_text())
    if manifest.get('schema')!='ttf_genetic_decker_pilot_geometry_v0.2':
        raise RuntimeError('canonical v0.2 geometry required')
    if ledger.get('schema')!='ttf_genetic_decker_pilot_geometry_source_v0.2':
        raise RuntimeError('canonical v0.2 source ledger required')
    if successor.get('schema')!='ttf_genetic_ibd_successor_v03_rule':
        raise RuntimeError('v0.3 successor rule required')
    if auth.get('schema')!='ttf_genetic_decker_execution_benchmark_authorization_v0.3':
        raise RuntimeError('v0.3 authorization drift')
    if auth.get('status')!='authorize_exactly_one_v03_successor_execution_benchmark':
        raise RuntimeError('v0.3 benchmark not authorized')
    if ledger['geometry_csv_sha256']!=sha256_path(args.geometry):
        raise RuntimeError('geometry CSV drift')
    if ledger['geometry_fingerprint_sha256']!=manifest['geometry_fingerprint_sha256']:
        raise RuntimeError('geometry fingerprint ledger drift')
    if ledger['rule_sha256']!=sha256_path(args.rule):
        raise RuntimeError('pilot rule drift')
    if auth['successor_rule_blob_sha']!=auth.get('frozen_successor_rule_blob_sha'):
        raise RuntimeError('successor rule authorization mismatch')

    geometries=load_geometries(args.geometry)
    fingerprint=geometry_fingerprint([
        SpeciesGeometry(species=n,coordinates=geometries[n].coordinates)
        for n in sorted(geometries)
    ])
    if fingerprint!=ledger['geometry_fingerprint_sha256']:
        raise RuntimeError('reconstructed geometry fingerprint drift')

    split=manifest['split']
    t0=time.perf_counter()
    design=prepare_genetic_ttf_design(
        geometries,
        train_species=split['train_species'],
        eval_species=split['eval_species'],
        bandwidth=float(rule['core_v011_inheritance']['bandwidth_km']),
        prior_strength=float(rule['core_v011_inheritance']['prior_strength']),
        segment_points=int(rule['core_v011_inheritance']['segment_points']),
        min_training_edges=int(rule['geometry_filter']['min_endpoint_disjoint_ibd_training_edges']),
        strength_neighbours=int(rule['core_v011_inheritance']['strength_neighbours']),
    )
    design_seconds=time.perf_counter()-t0

    case_results=[]
    for case in auth['authorized_cases']:
        world=simulate_genetic_distance_world(
            geometries,
            shared_fraction=float(case['shared_fraction']),
            residual_amplitude=float(case['residual_amplitude']),
            ibd_strength=float(rule['genetic_world']['ibd_strength_primary']),
            noise_sd=float(case['noise_sd']),
            transition_width=float(rule['genetic_world']['transition_width']),
            noise_dimensions=int(rule['genetic_world']['noise_dimensions']),
            seed=int(case['seed']),
        )
        started=time.perf_counter()
        score=score_genetic_world(design,world)
        elapsed=time.perf_counter()-started
        finite=sum(np.isfinite(v) for v in score.species_scores.values())
        case_results.append({
            'label':str(case['label']),
            'shared_fraction':float(case['shared_fraction']),
            'residual_amplitude':float(case['residual_amplitude']),
            'noise_sd':float(case['noise_sd']),
            'statistic':float(score.statistic),
            'training_strength':float(score.training_strength),
            'finite_eval_species_scores':int(finite),
            'score_seconds':float(elapsed),
        })

    inv=auth['required_invariants']
    ibd=next(row for row in case_results if row['label']=='ibd_only')
    if abs(ibd['statistic'])>float(inv['ibd_only_statistic_absolute_max']):
        raise RuntimeError('v0.3 IBD-only heldout statistic invariant failed')
    if abs(ibd['training_strength'])>float(inv['ibd_only_training_strength_absolute_max']):
        raise RuntimeError('v0.3 IBD-only training-strength invariant failed')
    if bool(inv['all_110_eval_species_scores_finite']) and any(
        row['finite_eval_species_scores']!=110 for row in case_results
    ):
        raise RuntimeError('v0.3 finite-score invariant failed')

    out={
        'schema':'ttf_genetic_decker_execution_benchmark_v0.3',
        'status':'successor_execution_benchmark_only_not_qualification',
        'parent_failure':'benchmarks/frozen/genetic_decker_execution_benchmark_v02_failure.json',
        'diagnosis':'benchmarks/frozen/genetic_ibd_only_numeric_diagnostic_v0.1.json',
        'successor_rule':str(args.successor_rule),
        'geometry_fingerprint_sha256':fingerprint,
        'species_count':len(geometries),
        'train_species':len(design.train_species),
        'eval_species':len(design.eval_species),
        'total_edges':int(sum(e.n_edges for e in design.template_edges.values())),
        'design_preparation_seconds':float(design_seconds),
        'cases':case_results,
        'empirical_genetic_outcomes_opened':False,
        'qualification_claim_made':False,
        'v02_failure_remains_immutable':True,
    }
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(out,indent=2,sort_keys=True)+'\n')
    print(json.dumps(out,sort_keys=True))
    return 0

if __name__=='__main__':
    raise SystemExit(main())
