#!/usr/bin/env python3
"""Post-exposure exploratory LGM design census; never reads genetic responses.

The rectangle is the original pilot's unsnapped coordinate bbox + 0.12 deg,
inside its outward-snapped raster window. Exact recorded-coordinate equality
is NOT evidence of equal coordinate precision, collection date or microhabitat.
Species labels in a pair are alphabetical, not independent source/target roles.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from collections import defaultdict
from decimal import Decimal, InvalidOperation
from itertools import combinations
from pathlib import Path

import numpy as np
from scipy.stats import rankdata

SOURCES = {
    'historical_host_memory_candidate_642_v0.1.csv':
        'c36cbb2ba0cbf7d0222645a04538c78236cfda392dd0a3d11dd443f12347d35b',
    'historical_host_memory_localities_v0.1.csv':
        '037cd8fa1f059fb67c349a465540d3d5fac469b5d14a2a4d2658d74c036c0ae9',
    'historical_host_climate_lcp_bungalotis_64edges_v0.1.csv':
        '7d255b208a485e7bc6d02e4e5748781c5e10170ed8bf161542dce3135c897b3b',
}
ANCHOR = 'Bungalotis erythus'


def verify_sha(path: Path, expected: str) -> None:
    if hashlib.sha256(path.read_bytes()).hexdigest() != expected:
        raise ValueError(f'SHA mismatch: {path.name}')


def read_csv(path: Path) -> list[dict]:
    with path.open(encoding='utf-8', newline='') as stream:
        return list(csv.DictReader(stream))


def host_ids(value: str) -> set[str]:
    result = set()
    for item in value.split(';'):
        try:
            number = Decimal(item.strip())
        except InvalidOperation as exc:
            raise ValueError('invalid host identifier') from exc
        if not number.is_finite() or number <= 0 or number != number.to_integral_value():
            raise ValueError('host identifiers must be finite positive integers')
        result.add(str(int(number)))
    return result


def coordinate(row: dict) -> tuple[Decimal, Decimal]:
    try:
        lat, lon = Decimal(str(row['latitude'])), Decimal(str(row['longitude']))
    except InvalidOperation as exc:
        raise ValueError('invalid coordinate') from exc
    if not lat.is_finite() or not lon.is_finite() or abs(lat) > 90 or abs(lon) > 180:
        raise ValueError('invalid coordinate')
    return lat.normalize(), lon.normalize()


def census(candidates: list[dict], localities: list[dict], bounds: tuple,
           minimum_localities: int = 12) -> dict:
    """Report all candidates; subset only on recorded geometry, never responses."""
    west, south, east, north = map(lambda x: Decimal(str(x)), bounds)
    if not (-180 <= west < east <= 180 and -90 <= south < north <= 90):
        raise ValueError('invalid non-wrapping census bounds')
    if minimum_localities < 1:
        raise ValueError('minimum_localities must be positive')
    by_species = {r['species']: r for r in candidates}
    if len(by_species) != len(candidates):
        raise ValueError('duplicate candidate species')
    points, inside = defaultdict(dict), defaultdict(dict)
    indices = defaultdict(set)
    for row in localities:
        sp = row['species']
        if sp not in by_species:
            raise ValueError('locality species missing from candidates')
        point = coordinate(row)
        idx = int(row['locality_index'])
        if idx < 0 or idx in indices[sp]:
            raise ValueError('duplicate or negative locality index')
        if point in points[sp]:
            raise ValueError('duplicate coordinate within species')
        points[sp][point] = idx
        indices[sp].add(idx)
        lat, lon = point
        if south <= lat <= north and west <= lon <= east:
            inside[sp][point] = idx
    screen, hosts = [], {}
    for sp, row in sorted(by_species.items()):
        hosts[sp] = host_ids(row['accepted_host_ids'])
        if len(hosts[sp]) != int(row['host_species']):
            raise ValueError('candidate host count mismatch')
        if 'n_localities' in row and len(points[sp]) != int(row['n_localities']):
            raise ValueError('candidate locality count mismatch')
        screen.append({
            'species': sp, 'family_as_frozen': row['family'],
            'localities_total': len(points[sp]), 'localities_inside': len(inside[sp]),
            'eligible_geometry': len(inside[sp]) >= minimum_localities,
            'host_species': len(hosts[sp]), 'host_families': int(row['host_families']),
            'accepted_host_ids': ';'.join(sorted(hosts[sp], key=int)),
            'accepted_host_names': row['accepted_host_names'].strip(),
        })
    panel = [r for r in screen if r['eligible_geometry']]
    pairs = []
    shared_index = {}
    for a, b in combinations([r['species'] for r in panel], 2):
        common = sorted(inside[a].keys() & inside[b].keys())
        paired = [{
            'latitude': str(lat), 'longitude': str(lon),
            'source_locality_index': inside[a][(lat, lon)],
            'target_locality_index': inside[b][(lat, lon)],
        } for lat, lon in common]
        # Decimal string pairs are stable identifiers without rounding/tolerance.
        key = a + '|' + b
        shared_index[key] = paired
        pairs.append({
            'source_species': a, 'target_species': b,
            'shared_recorded_coordinates': len(common),
            'same_frozen_family': by_species[a]['family'] == by_species[b]['family'],
            'host_shared_species': len(hosts[a] & hosts[b]),
            'host_union_species': len(hosts[a] | hosts[b]),
            'host_jaccard': len(hosts[a] & hosts[b]) / len(hosts[a] | hosts[b]),
            'shared_coordinate_pair_queries': len(common) * (len(common)-1) // 2,
        })
    return {'screen': screen, 'panel': panel, 'pairs': pairs, 'shared_index': shared_index,
            'inside_species_site_records': sum(len(inside[r['species']]) for r in panel),
            'inside_distinct_recorded_coordinates': len(set().union(
                *(inside[r['species']].keys() for r in panel))) if panel else 0}


def shared_queries(result: dict, minimum_shared: int = 12,
                   identical_host_only: bool = False) -> list[dict]:
    if minimum_shared < 2:
        raise ValueError('at least two shared coordinates required')
    queries = []
    for pair in result['pairs']:
        if pair['shared_recorded_coordinates'] < minimum_shared:
            continue
        if identical_host_only and pair['host_jaccard'] != 1.0:
            continue
        a, b = pair['source_species'], pair['target_species']
        for left, right in combinations(result['shared_index'][a+'|'+b], 2):
            queries.append({
                'source_species': a, 'target_species': b,
                'latitude_left': left['latitude'], 'longitude_left': left['longitude'],
                'latitude_right': right['latitude'], 'longitude_right': right['longitude'],
                'source_node_left': left['source_locality_index'],
                'source_node_right': right['source_locality_index'],
                'target_node_left': left['target_locality_index'],
                'target_node_right': right['target_locality_index'],
            })
    return queries


def rank_correlation(x, y) -> float | None:
    a, b = np.asarray(x, dtype=float), np.asarray(y, dtype=float)
    if a.ndim != 1 or b.shape != a.shape or not np.isfinite(a).all() or not np.isfinite(b).all():
        raise ValueError('invalid correlation vectors')
    if len(a) < 3 or np.ptp(a) == 0 or np.ptp(b) == 0:
        return None
    return float(np.corrcoef(rankdata(a), rankdata(b))[0, 1])


def shape_diagnostic(rows: list[dict]) -> dict:
    rows = sorted(rows, key=lambda r: int(r['edge_index']))
    fields = ('R_host_0BP_km', 'R_host_LGM_km', 'R_self_0BP_km', 'R_self_LGM_km')
    costs = np.array([[float(r[f]) for f in fields] for r in rows])
    distance = np.array([float(r['geographic_chord_km']) for r in rows])
    if not rows or len({int(r['edge_index']) for r in rows}) != len(rows):
        raise ValueError('empty or duplicate edge data')
    if not np.isfinite(costs).all() or not np.isfinite(distance).all() or np.any(costs < 0) or np.any(distance <= 0):
        raise ValueError('invalid path cost or geographic distance')
    zero = (costs == 0).all(axis=1)
    if np.any((costs == 0).any(axis=1) & ~zero):
        raise ValueError('only some scenarios have a zero path')
    positive = ~zero
    delta_h = costs[:, 1] - costs[:, 0]
    delta_s = costs[:, 3] - costs[:, 2]
    h, s, g = delta_h[positive], delta_s[positive], distance[positive]
    return {
        'all_edges': len(rows), 'positive_path_edges': int(positive.sum()),
        'same_pixel_edge_ids': [int(r['edge_index']) for r, z in zip(rows, zero) if z],
        'mean_host_delta_all_edges_cost_km': float(delta_h.mean()),
        'mean_self_delta_all_edges_cost_km': float(delta_s.mean()),
        'spearman_delta_raw': rank_correlation(h, s),
        'spearman_delta_per_chord_km': rank_correlation(h/g, s/g),
        'median_host_LGM_to_0BP_ratio': float(np.median(costs[positive, 1]/costs[positive, 0])) if positive.any() else None,
        'median_self_LGM_to_0BP_ratio': float(np.median(costs[positive, 3]/costs[positive, 2])) if positive.any() else None,
        'p_values_computed': False,
        'interpretation': 'Amplitude contrast is not spatial-rank divergence. Normalizing by distance is a descriptive sensitivity, not complete geographic confound control.',
        'reference_counts': {'host': 665, 'self': 26},
        'equal_reference_count_edgewise_shape_checked': False,
    }


def write_csv(path: Path, rows: list[dict]) -> None:
    if not rows:
        raise ValueError('refusing an ambiguous headerless empty table')
    with path.open('w', newline='', encoding='utf-8') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]), lineterminator='\n')
        writer.writeheader()
        writer.writerows(rows)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--sources-directory', type=Path, required=True)
    parser.add_argument('--output-directory', type=Path, required=True)
    args = parser.parse_args()
    for name, sha in SOURCES.items():
        verify_sha(args.sources_directory/name, sha)
    candidates = read_csv(args.sources_directory/list(SOURCES)[0])
    localities = read_csv(args.sources_directory/list(SOURCES)[1])
    pilot = read_csv(args.sources_directory/list(SOURCES)[2])
    anchor = [coordinate(r) for r in localities if r['species'] == ANCHOR]
    padding = Decimal('0.12')
    bounds = (min(p[1] for p in anchor)-padding, min(p[0] for p in anchor)-padding,
              max(p[1] for p in anchor)+padding, max(p[0] for p in anchor)+padding)
    result = census(candidates, localities, bounds, minimum_localities=12)
    queries = shared_queries(result, minimum_shared=12)
    identity_queries = shared_queries(result, minimum_shared=2, identical_host_only=True)
    out = args.output_directory
    out.mkdir(parents=True, exist_ok=True)
    tables = {
        'all_642_species_geometry_screen.csv': result['screen'],
        'regional_species_panel.csv': result['panel'],
        'all_regional_species_pairs.csv': result['pairs'],
        'shared12_coordinate_pair_queries.csv': queries,
        'identical_host_coordinate_pair_queries.csv': identity_queries,
    }
    for name, table in tables.items():
        write_csv(out/name, table)
    pairs12 = [r for r in result['pairs'] if r['shared_recorded_coordinates'] >= 12]
    summary = {
        'schema': 'ttf_exploratory_same_site_host_landscape_design_v0.1',
        'status': 'POST_EXPOSURE_EXPLORATORY_DESIGN_CENSUS_NOT_LGM_EFFECT_TEST',
        'date': '2026-10-19', 'source_commit': '7319fc9b555aabb2bd6e543dec90ee632588e33c',
        'source_sha256': SOURCES,
        'census_bounds_west_south_east_north': [str(x) for x in bounds],
        'bounds_definition': 'Original anchor coordinate rectangle plus 0.12 degrees; BEFORE outward native-pixel snapping. Not the exact raster-cell extent.',
        'minimum_localities_in_rectangle': 12,
        'candidate_species': len(candidates), 'input_locality_rows': len(localities),
        'regional_species': len(result['panel']),
        'regional_families_as_frozen': len({r['family_as_frozen'] for r in result['panel']}),
        'regional_species_site_records': result['inside_species_site_records'],
        'regional_unique_recorded_coordinates': result['inside_distinct_recorded_coordinates'],
        'regional_unordered_species_pairs': len(result['pairs']),
        'shared_coordinate_threshold_diagnostics': {
            str(n): {'species_pairs': sum(p['shared_recorded_coordinates'] >= n for p in result['pairs']),
                     'host_disjoint_pairs': sum(p['shared_recorded_coordinates'] >= n and p['host_shared_species']==0 for p in result['pairs'])}
            for n in (1, 3, 6, 12)},
        'shared12_pairs': pairs12,
        'shared12_unique_species': len({p[k] for p in pairs12 for k in ('source_species','target_species')}),
        'shared12_pair_query_rows': len(queries),
        'identical_host_pairs': [p for p in result['pairs'] if p['host_jaccard']==1.0],
        'identical_host_query_rows': len(identity_queries),
        'pilot_shape_diagnostic': shape_diagnostic(pilot),
        'output_sha256': {name: hashlib.sha256((out/name).read_bytes()).hexdigest() for name in tables},
        'boundaries': {
            'new_nucleotide_identity_read': False, 'genetic_test_run': False,
            'new_multispecies_climate_routing_run': False,
            'native_raster_validity_and_endpoint_distinctness_pending': True,
            'shared_recorded_coordinates_prove_same_microhabitat': False,
            'query_pairs_independent_replicates': False,
            'species_disjoint_or_source_disjoint_confirmation': False,
            'original_259_state': 'NOT_EVALUABLE_EMPIRICAL_RESPONSE',
            'regional_subset_selected_after_previous_climate_exposure': True,
            'complete_host_knowledge_assumed': False,
            'all_642_retained_in_screen_with_exclusions_visible': True,
        },
    }
    (out/'same_site_design_summary.json').write_text(json.dumps(summary, indent=2, ensure_ascii=False, allow_nan=False)+'\n', encoding='utf-8')
    print(json.dumps({k:summary[k] for k in ('status','regional_species','regional_unordered_species_pairs','shared12_unique_species','shared12_pair_query_rows')}, indent=2))
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
