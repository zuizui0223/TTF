#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
from pathlib import Path

from benchmark_genetic_decker_pipeline import load_geometries
from freeze_decker_genetic_pilot_geometry import nearest_training_support, species_edges


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument('--geometry', type=Path, required=True)
    ap.add_argument('--manifest', type=Path, required=True)
    ap.add_argument('--pilot-rule', type=Path, required=True)
    ap.add_argument('--output', type=Path, required=True)
    args = ap.parse_args()

    manifest = json.loads(args.manifest.read_text())
    rule = json.loads(args.pilot_rule.read_text())
    if manifest.get('schema') != 'ttf_genetic_decker_pilot_geometry_v0.2':
        raise RuntimeError('canonical frozen geometry required')
    if manifest.get('genetic_outcomes_opened') is not False:
        raise RuntimeError('empirical outcome firewall is open')
    if any(bool(v) for v in rule['outcome_firewall'].values()):
        raise RuntimeError('pilot-rule outcome firewall is open')

    geometries = load_geometries(args.geometry)
    panels = manifest['selected_panels']
    classes = {name: str(panels[name]['class']) for name in geometries}
    aves = tuple(sorted(name for name, cls in classes.items() if cls == 'Aves'))
    bats = tuple(sorted(name for name, cls in classes.items() if cls == 'Chiroptera'))
    if len(aves) != 149 or len(bats) != 72:
        raise RuntimeError(f'class-count drift: Aves={len(aves)}, Chiroptera={len(bats)}')

    edges = {
        name: species_edges(
            name,
            geometries[name].coordinates,
            geometries[name].edge_nodes,
        )
        for name in geometries
    }
    bandwidth = float(rule['core_v011_inheritance']['bandwidth_km'])

    def direction(train: tuple[str, ...], evaluation: tuple[str, ...]) -> dict:
        result = nearest_training_support(
            {name: edges[name] for name in train},
            {name: edges[name] for name in evaluation},
            bandwidth_km=bandwidth,
        )
        return result

    out = {
        'schema': 'ttf_genetic_reciprocal_class_support_census_v0.1',
        'status': 'response_blind_geometry_only',
        'geometry_fingerprint_sha256': manifest['geometry_fingerprint_sha256'],
        'bandwidth_km': bandwidth,
        'class_counts': {'Aves': len(aves), 'Chiroptera': len(bats)},
        'directions': {
            'Aves_to_Chiroptera': direction(aves, bats),
            'Chiroptera_to_Aves': direction(bats, aves),
        },
        'empirical_genetic_outcomes_opened': False,
        'qualification_claim_made': False,
        'claim_boundary': 'This census evaluates geometry/support only. Reciprocal class transfer is not inferentially authorized until each direction receives its own synthetic qualification on the exact class split.'
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(out, indent=2, sort_keys=True) + '\n')
    compact = {
        'schema': out['schema'],
        'class_counts': out['class_counts'],
        'bandwidth_km': bandwidth,
        'directions': {
            key: {k:v for k,v in value.items() if k != 'per_eval_species'}
            for key, value in out['directions'].items()
        },
        'empirical_genetic_outcomes_opened': False,
    }
    print(json.dumps(compact, sort_keys=True))
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
