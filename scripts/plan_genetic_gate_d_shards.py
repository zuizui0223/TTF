#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
from pathlib import Path


def make_ranges(total: int, shard_size: int) -> list[tuple[int,int]]:
    if total < 1 or shard_size < 1:
        raise ValueError('total and shard size must be positive')
    return [(start, min(start + shard_size, total)) for start in range(0, total, shard_size)]


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument('--pilot-rule', type=Path, required=True)
    ap.add_argument('--authorization', type=Path, required=True)
    ap.add_argument('--output', type=Path, required=True)
    args = ap.parse_args()

    pilot = json.loads(args.pilot_rule.read_text())
    auth = json.loads(args.authorization.read_text())
    if auth.get('schema') != 'ttf_genetic_gate_d_execution_authorization_v0.2':
        raise RuntimeError('Gate-D execution authorization drift')
    ref_total = int(pilot['qualification']['reference_worlds_per_configuration'])
    obs_total = int(pilot['qualification']['observed_worlds_per_cell'])
    ref_size = int(auth['reference_shard_size'])
    obs_size = int(auth['observed_shard_size'])

    reference_include = []
    for label in pilot['genetic_world']['private_reference_configurations']:
        for start, stop in make_ranges(ref_total, ref_size):
            reference_include.append({
                'config': str(label),
                'start': int(start),
                'stop': int(stop),
                'artifact': f'genetic-gate-ref-{label}-{start}-{stop}',
            })
    observed_include = []
    for shared, amplitude in pilot['qualification']['mandatory_primary_cells']:
        shared = float(shared)
        amplitude = float(amplitude)
        cell = f's{shared:g}-a{amplitude:g}'.replace('.', 'p')
        for start, stop in make_ranges(obs_total, obs_size):
            observed_include.append({
                'shared': shared,
                'amplitude': amplitude,
                'start': int(start),
                'stop': int(stop),
                'cell': cell,
                'artifact': f'genetic-gate-obs-{cell}-{start}-{stop}',
            })

    out = {
        'schema': 'ttf_genetic_gate_d_shard_plan_v0.1',
        'status': 'deterministic_execution_plan',
        'reference_matrix': {'include': reference_include},
        'observed_matrix': {'include': observed_include},
        'reference_job_count': len(reference_include),
        'observed_job_count': len(observed_include),
        'reference_worlds_per_configuration': ref_total,
        'observed_worlds_per_cell': obs_total,
        'selected_batch_width': int(auth['selected_batch_width']),
        'empirical_genetic_outcomes_opened': False,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(out, indent=2, sort_keys=True) + '\n')
    print(json.dumps(out, sort_keys=True))
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
