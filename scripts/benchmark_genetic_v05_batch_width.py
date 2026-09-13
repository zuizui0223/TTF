#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import time
from pathlib import Path

import numpy as np

from benchmark_genetic_decker_pipeline import load_geometries
from ttf.calibration import seed_for
from ttf.genetic_batch_execution import (
    prepare_genetic_cached_transfer,
    score_genetic_world_batch,
)
from ttf.genetic_gate import prepare_genetic_ttf_design, score_genetic_world
from ttf.genetic_simulate import simulate_genetic_distance_world


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument('--geometry', type=Path, required=True)
    ap.add_argument('--manifest', type=Path, required=True)
    ap.add_argument('--pilot-rule', type=Path, required=True)
    ap.add_argument('--authorization', type=Path, required=True)
    ap.add_argument('--batch-width', type=int, required=True)
    ap.add_argument('--output', type=Path, required=True)
    args = ap.parse_args()

    manifest = json.loads(args.manifest.read_text())
    pilot = json.loads(args.pilot_rule.read_text())
    auth = json.loads(args.authorization.read_text())
    if auth.get('schema') != 'ttf_genetic_v05_batch_width_benchmark_authorization':
        raise RuntimeError('v0.5 batch benchmark authorization drift')
    if auth.get('status') != 'authorize_runtime_only_batch_width_benchmark':
        raise RuntimeError('v0.5 batch benchmark not authorized')
    width = int(args.batch_width)
    if width not in [int(v) for v in auth['candidate_batch_widths']]:
        raise RuntimeError('unauthorized batch width')

    geometries = load_geometries(args.geometry)
    split = manifest['split']
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

    family = auth['world_family']
    n_worlds = int(auth['n_worlds'])
    worlds = [
        simulate_genetic_distance_world(
            geometries,
            shared_fraction=float(family['shared_fraction']),
            residual_amplitude=float(family['residual_amplitude']),
            ibd_strength=float(pilot['genetic_world']['ibd_strength_primary']),
            noise_sd=float(family['noise_sd']),
            transition_width=float(pilot['genetic_world']['transition_width']),
            noise_dimensions=int(pilot['genetic_world']['noise_dimensions']),
            seed=seed_for(
                int(family['master_seed']),
                'genetic_v05_batch_benchmark',
                replicate,
            ),
        )
        for replicate in range(n_worlds)
    ]

    sequential = score_genetic_world(design, worlds[0])

    cache_started = time.perf_counter()
    cached = prepare_genetic_cached_transfer(
        design,
        edge_chunk_size=32,
        train_chunk_size=4096,
    )
    cache_seconds = time.perf_counter() - cache_started

    statistics: list[float] = []
    strengths: list[float] = []
    species_scores = {name: [] for name in design.eval_species}
    n_eval: list[int] = []
    score_started = time.perf_counter()
    for start in range(0, n_worlds, width):
        stop = min(start + width, n_worlds)
        scored = score_genetic_world_batch(design, worlds[start:stop], cached)
        statistics.extend(map(float, scored.statistics))
        strengths.extend(map(float, scored.training_strengths))
        n_eval.extend(map(int, scored.n_eval_species))
        for name in design.eval_species:
            species_scores[name].extend(map(float, scored.species_scores[name]))
    score_seconds = time.perf_counter() - score_started

    if len(statistics) != n_worlds or len(strengths) != n_worlds or len(n_eval) != n_worlds:
        raise RuntimeError('batch output count drift')
    finite_invariant = (
        all(value == len(design.eval_species) for value in n_eval)
        and np.isfinite(np.asarray(statistics)).all()
        and np.isfinite(np.asarray(strengths)).all()
        and all(np.isfinite(np.asarray(values)).all() for values in species_scores.values())
    )
    if not finite_invariant:
        raise RuntimeError('finite-score invariant failed')

    statistic_difference = abs(float(statistics[0]) - float(sequential.statistic))
    strength_difference = abs(float(strengths[0]) - float(sequential.training_strength))
    species_difference = max(
        abs(float(species_scores[name][0]) - float(sequential.species_scores[name]))
        for name in design.eval_species
    )
    if max(statistic_difference, strength_difference, species_difference) > 1e-12:
        raise RuntimeError('real-geometry sequential equivalence failed')

    out = {
        'schema': 'ttf_genetic_v05_batch_width_benchmark_v0.1',
        'status': 'runtime_only_candidate_passed_equivalence',
        'batch_width': width,
        'n_worlds': n_worlds,
        'cache_preparation_seconds': float(cache_seconds),
        'total_scoring_seconds': float(score_seconds),
        'seconds_per_world': float(score_seconds / n_worlds),
        'max_abs_statistic_difference_from_sequential_first_world': float(statistic_difference),
        'max_abs_training_strength_difference_from_sequential_first_world': float(strength_difference),
        'max_abs_species_score_difference_from_sequential_first_world': float(species_difference),
        'finite_score_invariant': True,
        'empirical_genetic_outcomes_opened': False,
        'qualification_claim_made': False,
        'statistic_values_persisted': False,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(out, indent=2, sort_keys=True) + '\n')
    print(json.dumps(out, sort_keys=True))
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
