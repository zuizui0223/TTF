#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
from pathlib import Path


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument('--input-dir', type=Path, required=True)
    ap.add_argument('--pilot-rule', type=Path, required=True)
    ap.add_argument('--output', type=Path, required=True)
    args = ap.parse_args()

    pilot = json.loads(args.pilot_rule.read_text())
    expected_labels = tuple(pilot['genetic_world']['private_reference_configurations'].keys())
    expected_n = int(pilot['qualification']['reference_worlds_per_configuration'])
    grouped: dict[str, dict[int, tuple[float, float]]] = {label: {} for label in expected_labels}

    files = sorted(args.input_dir.rglob('*.json'))
    if not files:
        raise RuntimeError('no reference shard files found')
    for path in files:
        shard = json.loads(path.read_text())
        if shard.get('schema') != 'ttf_genetic_gate_d_reference_shard_v0.1':
            continue
        if shard.get('empirical_genetic_outcomes_opened') is not False:
            raise RuntimeError(f'empirical outcome firewall drift in {path}')
        label = str(shard['configuration'])
        if label not in grouped:
            raise RuntimeError(f'unexpected reference configuration {label}')
        for row in shard['rows']:
            replicate = int(row['replicate'])
            if replicate in grouped[label]:
                raise RuntimeError(f'duplicate reference replicate {label}:{replicate}')
            grouped[label][replicate] = (
                float(row['training_strength']),
                float(row['statistic']),
            )

    references: dict[str, dict[str, list[float]]] = {}
    for label in expected_labels:
        rows = grouped[label]
        if sorted(rows) != list(range(expected_n)):
            missing = sorted(set(range(expected_n)) - set(rows))[:10]
            raise RuntimeError(f'incomplete reference configuration {label}; first missing={missing}')
        references[label] = {
            'training_strength': [rows[i][0] for i in range(expected_n)],
            'statistic': [rows[i][1] for i in range(expected_n)],
        }

    out = {
        'schema': 'ttf_genetic_gate_d_references_v0.1',
        'status': 'ordered_private_reference_families_complete',
        'reference_worlds_per_configuration': expected_n,
        'profile_draws': int(pilot['qualification']['profile_strength_draws']),
        'calibration_statistic_draws': int(pilot['qualification']['calibration_statistic_draws']),
        'configuration_order': list(expected_labels),
        'references': references,
        'empirical_genetic_outcomes_opened': False,
        'qualification_claim_made': False,
    }
    if out['profile_draws'] + out['calibration_statistic_draws'] != expected_n:
        raise RuntimeError('profile/calibration draw split does not cover frozen references exactly')
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(out, indent=2, sort_keys=True) + '\n')
    print(json.dumps({k:v for k,v in out.items() if k != 'references'}, sort_keys=True))
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
