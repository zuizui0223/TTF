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
    ap = argparse.ArgumentParser()
    ap.add_argument('--geometry', type=Path, required=True)
    ap.add_argument('--manifest', type=Path, required=True)
    ap.add_argument('--source-ledger', type=Path, required=True)
    ap.add_argument('--pilot-rule', type=Path, required=True)
    ap.add_argument('--exact-fast-rule', type=Path, required=True)
    ap.add_argument('--exact-scale-addendum', type=Path, required=True)
    ap.add_argument('--authorization', type=Path, required=True)
    ap.add_argument('--output', type=Path, required=True)
    args = ap.parse_args()

    manifest = json.loads(args.manifest.read_text())
    ledger = json.loads(args.source_ledger.read_text())
    pilot = json.loads(args.pilot_rule.read_text())
    fast_rule = json.loads(args.exact_fast_rule.read_text())
    scale_rule = json.loads(args.exact_scale_addendum.read_text())
    auth = json.loads(args.authorization.read_text())

    if manifest.get('schema') != 'ttf_genetic_decker_pilot_geometry_v0.2':
        raise RuntimeError('canonical v0.2 geometry required')
    if ledger.get('schema') != 'ttf_genetic_decker_pilot_geometry_source_v0.2':
        raise RuntimeError('canonical v0.2 source ledger required')
    if fast_rule.get('schema') != 'ttf_genetic_ibd_exact_fast_v04_rule':
        raise RuntimeError('v0.4 exact-fast rule required')
    if scale_rule.get('schema') != 'ttf_genetic_ibd_v04_exact_scale_addendum':
        raise RuntimeError('v0.4 exact-scale addendum required')
    if auth.get('schema') != 'ttf_genetic_decker_execution_benchmark_authorization_v0.4':
        raise RuntimeError('v0.4 authorization drift')
    if auth.get('status') != 'authorize_exactly_one_v04_exact_fast_execution_benchmark':
        raise RuntimeError('v0.4 benchmark not authorized')
    if ledger['geometry_csv_sha256'] != sha256_path(args.geometry):
        raise RuntimeError('geometry CSV drift')
    if ledger['geometry_fingerprint_sha256'] != manifest['geometry_fingerprint_sha256']:
        raise RuntimeError('geometry fingerprint ledger drift')
    if ledger['rule_sha256'] != sha256_path(args.pilot_rule):
        raise RuntimeError('pilot rule drift')

    geometries = load_geometries(args.geometry)
    fingerprint = geometry_fingerprint([
        SpeciesGeometry(species=name, coordinates=geometries[name].coordinates)
        for name in sorted(geometries)
    ])
    if fingerprint != ledger['geometry_fingerprint_sha256']:
        raise RuntimeError('reconstructed geometry fingerprint drift')

    split = manifest['split']
    design_started = time.perf_counter()
    design = prepare_genetic_ttf_design(
        geometries,
        train_species=split['train_species'],
        eval_species=split['eval_species'],
        bandwidth=float(pilot['core_v011_inheritance']['bandwidth_km']),
        prior_strength=float(pilot['core_v011_inheritance']['prior_strength']),
        segment_points=int(pilot['core_v011_inheritance']['segment_points']),
        min_training_edges=int(pilot['geometry_filter']['min_endpoint_disjoint_ibd_training_edges']),
        strength_neighbours=int(pilot['core_v011_inheritance']['strength_neighbours']),
    )
    design_seconds = time.perf_counter() - design_started

    case_results = []
    total_started = time.perf_counter()
    for case in auth['authorized_cases']:
        world = simulate_genetic_distance_world(
            geometries,
            shared_fraction=float(case['shared_fraction']),
            residual_amplitude=float(case['residual_amplitude']),
            ibd_strength=float(pilot['genetic_world']['ibd_strength_primary']),
            noise_sd=float(case['noise_sd']),
            transition_width=float(pilot['genetic_world']['transition_width']),
            noise_dimensions=int(pilot['genetic_world']['noise_dimensions']),
            seed=int(case['seed']),
        )
        started = time.perf_counter()
        score = score_genetic_world(design, world)
        elapsed = time.perf_counter() - started
        finite = sum(np.isfinite(v) for v in score.species_scores.values())
        case_results.append({
            'label': str(case['label']),
            'shared_fraction': float(case['shared_fraction']),
            'residual_amplitude': float(case['residual_amplitude']),
            'noise_sd': float(case['noise_sd']),
            'statistic': float(score.statistic),
            'training_strength': float(score.training_strength),
            'finite_eval_species_scores': int(finite),
            'score_seconds': float(elapsed),
        })
    total_score_seconds = time.perf_counter() - total_started

    inv = auth['required_invariants']
    ibd = next(row for row in case_results if row['label'] == 'ibd_only')
    if abs(ibd['statistic']) > float(inv['ibd_only_statistic_absolute_max']):
        raise RuntimeError('v0.4 IBD-only heldout statistic invariant failed')
    if abs(ibd['training_strength']) > float(inv['ibd_only_training_strength_absolute_max']):
        raise RuntimeError('v0.4 IBD-only training-strength invariant failed')
    if bool(inv['all_110_eval_species_scores_finite']) and any(
        row['finite_eval_species_scores'] != 110 for row in case_results
    ):
        raise RuntimeError('v0.4 finite-score invariant failed')

    out = {
        'schema': 'ttf_genetic_decker_execution_benchmark_v0.4',
        'status': 'exact_fast_execution_benchmark_passed_invariants_not_qualification',
        'parent_failures': [
            'benchmarks/frozen/genetic_decker_execution_benchmark_v02_failure.json',
            'benchmarks/frozen/genetic_decker_execution_benchmark_v03_failure.json',
        ],
        'exact_fast_rule': str(args.exact_fast_rule),
        'exact_scale_addendum': str(args.exact_scale_addendum),
        'geometry_fingerprint_sha256': fingerprint,
        'species_count': len(geometries),
        'train_species': len(design.train_species),
        'eval_species': len(design.eval_species),
        'total_edges': int(sum(edge.n_edges for edge in design.template_edges.values())),
        'design_preparation_seconds': float(design_seconds),
        'total_score_seconds': float(total_score_seconds),
        'cases': case_results,
        'empirical_genetic_outcomes_opened': False,
        'qualification_claim_made': False,
        'v02_v03_failures_remain_immutable': True,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(out, indent=2, sort_keys=True) + '\n')
    print(json.dumps(out, sort_keys=True))
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
