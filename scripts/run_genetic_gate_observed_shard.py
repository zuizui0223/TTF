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
from ttf.profiled_private_null import profiled_private_pvalue


FIREWALL_KEYS = {
    'empirical_pairwise_genetic_distances_opened',
    'sequence_characters_opened_for_inference',
    'alignment_divergence_summaries_opened',
    'monmonier_results_opened',
    'break_presence_absence_opened',
}


def _assert_closed_firewall(auth: dict) -> None:
    firewall = auth.get('outcome_firewall')
    if not isinstance(firewall, dict) or set(firewall) != FIREWALL_KEYS:
        raise RuntimeError('empirical genetic outcome firewall schema drift')
    if any(value is not False for value in firewall.values()):
        raise RuntimeError('empirical genetic outcome firewall is open')


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument('--geometry', type=Path, required=True)
    ap.add_argument('--manifest', type=Path, required=True)
    ap.add_argument('--pilot-rule', type=Path, required=True)
    ap.add_argument('--authorization', type=Path, required=True)
    ap.add_argument('--references', type=Path, required=True)
    ap.add_argument('--shared-fraction', type=float, required=True)
    ap.add_argument('--residual-amplitude', type=float, required=True)
    ap.add_argument('--start', type=int, required=True)
    ap.add_argument('--stop', type=int, required=True)
    ap.add_argument('--output', type=Path, required=True)
    args = ap.parse_args()

    manifest = json.loads(args.manifest.read_text())
    pilot = json.loads(args.pilot_rule.read_text())
    auth = json.loads(args.authorization.read_text())
    refs_payload = json.loads(args.references.read_text())
    if auth.get('schema') != 'ttf_genetic_gate_d_execution_authorization_v0.2':
        raise RuntimeError('Gate-D execution authorization drift')
    if auth.get('status') != 'authorize_frozen_genetic_gate_d_synthetic_qualification':
        raise RuntimeError('Gate-D execution is not authorized')
    if refs_payload.get('schema') != 'ttf_genetic_gate_d_references_v0.1':
        raise RuntimeError('ordered reference payload required')
    _assert_closed_firewall(auth)

    shared = float(args.shared_fraction)
    amplitude = float(args.residual_amplitude)
    allowed_cells = {
        (float(cell[0]), float(cell[1]))
        for cell in pilot['qualification']['mandatory_primary_cells']
    }
    if (shared, amplitude) not in allowed_cells:
        raise RuntimeError('observed cell is not one of the frozen mandatory cells')
    total = int(pilot['qualification']['observed_worlds_per_cell'])
    start, stop = int(args.start), int(args.stop)
    if not 0 <= start < stop <= total:
        raise RuntimeError('observed shard range is outside frozen replicate range')
    batch_width = int(auth['selected_batch_width'])

    references = {
        label: (
            np.asarray(data['training_strength'], dtype=float),
            np.asarray(data['statistic'], dtype=float),
        )
        for label, data in refs_payload['references'].items()
    }
    expected_labels = set(pilot['genetic_world']['private_reference_configurations'])
    if set(references) != expected_labels:
        raise RuntimeError('reference configuration set drift')

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

    master_seed = int(auth['master_seed'])
    rows: list[dict] = []
    for batch_start in range(start, stop, batch_width):
        batch_stop = min(batch_start + batch_width, stop)
        worlds = [
            simulate_genetic_distance_world(
                geometries,
                shared_fraction=shared,
                residual_amplitude=amplitude,
                ibd_strength=float(pilot['genetic_world']['ibd_strength_primary']),
                noise_sd=float(pilot['genetic_world']['locality_noise_sd']),
                transition_width=float(pilot['genetic_world']['transition_width']),
                noise_dimensions=int(pilot['genetic_world']['noise_dimensions']),
                seed=seed_for(
                    master_seed,
                    'genetic_observed',
                    shared,
                    amplitude,
                    replicate,
                ),
            )
            for replicate in range(batch_start, batch_stop)
        ]
        scored = score_genetic_world_batch(design, worlds, cached)
        for column, replicate in enumerate(range(batch_start, batch_stop)):
            if int(scored.n_eval_species[column]) != len(design.eval_species):
                raise RuntimeError('non-finite held-out species count in observed shard')
            p_value, selected, component, distances = profiled_private_pvalue(
                float(scored.statistics[column]),
                float(scored.training_strengths[column]),
                references,
                profile_draws=int(pilot['qualification']['profile_strength_draws']),
                selected_configs=int(pilot['core_v011_inheritance']['profiled_private_selected_configurations']),
            )
            rows.append({
                'replicate': int(replicate),
                'p_value': float(p_value),
                'statistic': float(scored.statistics[column]),
                'training_strength': float(scored.training_strengths[column]),
                'selected_configurations': list(selected),
                'component_p_values': {k: float(v) for k,v in component.items()},
                'profile_distances': {k: float(v) for k,v in distances.items()},
            })

    if [row['replicate'] for row in rows] != list(range(start, stop)):
        raise RuntimeError('observed shard replicate ordering drift')
    out = {
        'schema': 'ttf_genetic_gate_d_observed_shard_v0.1',
        'status': 'synthetic_mandatory_observed_cell_shard',
        'shared_fraction': shared,
        'residual_amplitude': amplitude,
        'start': start,
        'stop': stop,
        'n_worlds': stop - start,
        'alpha': float(pilot['qualification']['alpha']),
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
