#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import time
from pathlib import Path

import numpy as np

from benchmark_genetic_decker_pipeline import load_geometries
from ttf.chunked_transfer import score_chunked_batch
from ttf.genetic_gate import prepare_genetic_ttf_design
from ttf.genetic_ibd import crossfit_ibd_residuals_prepared
from ttf.genetic_simulate import simulate_genetic_distance_world
from ttf.geometry_control import length_orthogonalized_turnover
from ttf.private_strength import training_private_strength_from_indices


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument('--geometry', type=Path, required=True)
    ap.add_argument('--manifest', type=Path, required=True)
    ap.add_argument('--pilot-rule', type=Path, required=True)
    ap.add_argument('--authorization', type=Path, required=True)
    ap.add_argument('--output', type=Path, required=True)
    args = ap.parse_args()

    manifest = json.loads(args.manifest.read_text())
    pilot = json.loads(args.pilot_rule.read_text())
    auth = json.loads(args.authorization.read_text())
    if auth.get('schema') != 'ttf_genetic_v04_execution_profile_authorization':
        raise RuntimeError('execution profile authorization drift')
    if auth.get('status') != 'authorize_one_computational_profile_after_v04_execution_pass':
        raise RuntimeError('execution profile not authorized')
    if auth.get('empirical_genetic_outcomes_opened') is not False:
        raise RuntimeError('empirical outcome firewall is open')

    geometries = load_geometries(args.geometry)
    split = manifest['split']
    started = time.perf_counter()
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
    design_seconds = time.perf_counter() - started

    case = auth['profile_world']
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

    used = design.train_species + design.eval_species
    train_response: dict[str, np.ndarray] = {}
    eval_response: dict[str, np.ndarray] = {}
    per_species = []
    ibd_total = 0.0
    control_total = 0.0

    total_started = time.perf_counter()
    for name in used:
        t0 = time.perf_counter()
        result = crossfit_ibd_residuals_prepared(
            world.genetic_distance[name],
            design.ibd_designs[name],
        )
        ibd_seconds = time.perf_counter() - t0
        ibd_total += ibd_seconds
        control_seconds = 0.0
        if name in design.train_species:
            t1 = time.perf_counter()
            train_response[name] = length_orthogonalized_turnover(
                result.residual_turnover,
                design.template_edges[name].length,
            )
            control_seconds = time.perf_counter() - t1
            control_total += control_seconds
        else:
            eval_response[name] = result.residual_turnover
        per_species.append({
            'species': name,
            'edges': int(design.template_edges[name].n_edges),
            'ibd_seconds': float(ibd_seconds),
            'training_control_seconds': float(control_seconds),
        })

    t2 = time.perf_counter()
    strength, _ = training_private_strength_from_indices(
        train_response,
        design.strength_indices,
        design.train_species,
    )
    strength_seconds = time.perf_counter() - t2

    t3 = time.perf_counter()
    scored = score_chunked_batch(
        design.prepared,
        {name: train_response[name][:, None] for name in design.train_species},
        {name: eval_response[name][:, None] for name in design.eval_species},
        edge_chunk_size=32,
        train_chunk_size=4096,
    )
    transfer_seconds = time.perf_counter() - t3
    total_seconds = time.perf_counter() - total_started

    per_species_sorted = sorted(
        per_species,
        key=lambda row: (-row['ibd_seconds'], row['species']),
    )
    out = {
        'schema': 'ttf_genetic_v04_execution_profile_v0.1',
        'status': 'computational_profile_only_not_qualification',
        'design_preparation_seconds': float(design_seconds),
        'profile_world': case,
        'stage_seconds': {
            'ibd_total': float(ibd_total),
            'training_geometry_control_total': float(control_total),
            'private_strength': float(strength_seconds),
            'transfer_score': float(transfer_seconds),
            'total_scoring_wall': float(total_seconds),
        },
        'stage_fraction_of_scoring_wall': {
            'ibd_total': float(ibd_total / total_seconds),
            'training_geometry_control_total': float(control_total / total_seconds),
            'private_strength': float(strength_seconds / total_seconds),
            'transfer_score': float(transfer_seconds / total_seconds),
        },
        'top_15_ibd_species': per_species_sorted[:15],
        'all_species_ibd_seconds_sum': float(sum(row['ibd_seconds'] for row in per_species)),
        'species_count': len(used),
        'total_edges': int(sum(design.template_edges[name].n_edges for name in used)),
        'finite_eval_species_scores': int(sum(np.isfinite(scored.species_scores[name][0]) for name in design.eval_species)),
        'training_strength_computed': bool(np.isfinite(strength)),
        'empirical_genetic_outcomes_opened': False,
        'qualification_claim_made': False,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(out, indent=2, sort_keys=True) + '\n')
    print(json.dumps(out, sort_keys=True))
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
