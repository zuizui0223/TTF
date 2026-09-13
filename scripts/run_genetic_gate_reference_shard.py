#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

from benchmark_genetic_decker_pipeline import load_geometries
from ttf.calibration import seed_for
from ttf.genetic_batch_execution import (
    prepare_genetic_cached_transfer,
    score_genetic_world_batch,
)
from ttf.genetic_gate import prepare_genetic_ttf_design
from ttf.genetic_simulate import simulate_genetic_distance_world


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument('--geometry', type=Path, required=True)
    ap.add_argument('--manifest', type=Path, required=True)
    ap.add_argument('--pilot-rule', type=Path, required=True)
    ap.add_argument('--authorization', type=Path, required=True)
    ap.add_argument('--config-label', required=True)
    ap.add_argument('--start', type=int, required=True)
    ap.add_argument('--stop', type=int, required=True)
    ap.add_argument('--output', type=Path, required=True)
    args = ap.parse_args()

    manifest = json.loads(args.manifest.read_text())
    pilot = json.loads(args.pilot_rule.read_text())
    auth = json.loads(args.authorization.read_text())
    if auth.get('schema') != 'ttf_genetic_gate_d_execution_authorization_v0.1':
        raise RuntimeError('Gate-D execution authorization drift')
    if auth.get('status') != 'authorize_frozen_genetic_gate_d_synthetic_qualification':
        raise RuntimeError('Gate-D execution is not authorized')
    if auth.get('empirical_genetic_outcomes_opened') is not False:
        raise RuntimeError('empirical genetic outcome firewall is open')

    label = str(args.config_label)
    configs = pilot['genetic_world']['private_reference_configurations']
    if label not in configs:
        raise RuntimeError(f'unknown private reference configuration: {label}')
    total = int(pilot['qualification']['reference_worlds_per_configuration'])
    start, stop = int(args.start), int(args.stop)
    if not 0 <= start < stop <= total:
        raise RuntimeError('reference shard range is outside frozen replicate range')
    batch_width = int(auth['selected_batch_width'])
    if batch_width < 1:
        raise RuntimeError('invalid authorized batch width')

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
    cached = prepare_genetic_cached_transfer(
        design,
        edge_chunk_size=int(auth['edge_chunk_size']),
        train_chunk_size=int(auth['train_chunk_size']),
    )

    amplitude, noise_sd = map(float, configs[label])
    master_seed = int(auth['master_seed'])
    rows: list[dict[str, float | int]] = []
    for batch_start in range(start, stop, batch_width):
        batch_stop = min(batch_start + batch_width, stop)
        worlds = [
            simulate_genetic_distance_world(
                geometries,
                shared_fraction=0.0,
                residual_amplitude=amplitude,
                ibd_strength=float(pilot['genetic_world']['ibd_strength_primary']),
                noise_sd=noise_sd,
                transition_width=float(pilot['genetic_world']['transition_width']),
                noise_dimensions=int(pilot['genetic_world']['noise_dimensions']),
                seed=seed_for(
                    master_seed,
                    'genetic_private_reference',
                    label,
                    replicate,
                ),
            )
            for replicate in range(batch_start, batch_stop)
        ]
        scored = score_genetic_world_batch(design, worlds, cached)
        for column, replicate in enumerate(range(batch_start, batch_stop)):
            if int(scored.n_eval_species[column]) != len(design.eval_species):
                raise RuntimeError('non-finite held-out species count in reference shard')
            rows.append({
                'replicate': int(replicate),
                'training_strength': float(scored.training_strengths[column]),
                'statistic': float(scored.statistics[column]),
            })

    if [row['replicate'] for row in rows] != list(range(start, stop)):
        raise RuntimeError('reference shard replicate ordering drift')
    out = {
        'schema': 'ttf_genetic_gate_d_reference_shard_v0.1',
        'status': 'synthetic_private_reference_shard',
        'configuration': label,
        'residual_amplitude': amplitude,
        'noise_sd': noise_sd,
        'start': start,
        'stop': stop,
        'n_worlds': stop - start,
        'rows': rows,
        'empirical_genetic_outcomes_opened': False,
        'qualification_claim_made': False,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(out, indent=2, sort_keys=True) + '\n')
    print(json.dumps({k:v for k,v in out.items() if k != 'rows'}, sort_keys=True))
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
