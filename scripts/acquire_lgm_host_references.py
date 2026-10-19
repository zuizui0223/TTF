#!/usr/bin/env python3
"""Acquire the frozen host taxon list once, without climate or genetic outcomes.

Raw quality-accepted observations are saved BEFORE spatial thinning. Spatial
exclusion must later be applied separately for each linked insect, using that
insect's exact evaluated geometry. This output is not a native-range, historical
occupancy, or paired-climate result. Successful and failed attempts are immutable.
"""
from __future__ import annotations
import argparse
import csv
import hashlib
import json
import os
import tempfile
import time
from collections import Counter
from concurrent.futures import ThreadPoolExecutor, wait, FIRST_COMPLETED
from datetime import datetime, timezone
from decimal import Decimal, InvalidOperation
from pathlib import Path
from typing import Any, Callable

import acquire_lgm_external_references as legacy

MANIFEST_SHA256 = 'a28fba9efa8a9dc79c162336f62c883059c5dddbfabec2861f554d038dd98044'
LEGACY_SCRIPT_SHA256 = '258553bf6b9180fdd03839a092fde2b890f0d77f2839f0ea7e5ff3c1ff64aaf2'
ALLOWED_RANKS = {'SPECIES', 'SUBSPECIES'}
STATUSES = {'OBSERVATIONS_ACQUIRED', 'TAXONOMIC_RANK_REVIEW', 'TAXONOMY_NOT_EXACT', 'REQUEST_ERROR'}


def normalize_id(value: Any) -> str:
    text = str(value).strip()
    try:
        number = Decimal(text)
    except InvalidOperation as exc:
        raise ValueError('invalid WCVP ID') from exc
    if (not number.is_finite() or number <= 0 or number != number.to_integral_value()
            or number.adjusted() > 12):
        raise ValueError('WCVP ID must be a positive bounded integer')
    return str(int(number))


def request_rank(name: str) -> str:
    tokens = name.split()
    if len(tokens) == 2:
        return 'SPECIES'
    if len(tokens) == 4 and tokens[2] == 'subsp.':
        return 'SUBSPECIES'
    return 'REQUIRES_TAXONOMIC_REVIEW'


def host_plan(rows: list[dict]) -> list[dict]:
    mapping: dict[str, dict] = {}
    for row in rows:
        key = normalize_id(row['wcvp_id_as_frozen'])
        name = row['accepted_host_name'].strip()
        insect = row['insect_species'].strip()
        if not name or not insect:
            raise ValueError('empty host/insect name')
        entry = mapping.setdefault(key, {'wcvp_id': key, 'name': name,
            'rank_hint': request_rank(name), 'insect_species': set()})
        if entry['name'] != name:
            raise ValueError('conflicting names for one frozen WCVP ID')
        entry['insect_species'].add(insect)
    result = [{**entry, 'insect_species': sorted(entry['insect_species']),
               'linked_insect_count': len(entry['insect_species'])} for entry in mapping.values()]
    return sorted(result, key=lambda entry: (-entry['linked_insect_count'], entry['name'], entry['wcvp_id']))


def plant_key(match: dict, name: str, rank: str) -> int | None:
    if (rank not in ALLOWED_RANKS or match.get('matchType') != 'EXACT'
            or match.get('kingdom') != 'Plantae' or match.get('rank') != rank
            or match.get('canonicalName', '').strip().casefold() != name.strip().casefold()
            or match.get('confidence', 0) < 95):
        return None
    value = match.get('acceptedUsageKey') or match.get('usageKey')
    return int(value) if value is not None else None


def clean_host_record(record: dict, key: int, rank: str) -> tuple[dict | None, str]:
    if rank == 'SPECIES':
        matched = record.get('speciesKey') == key
    elif rank == 'SUBSPECIES':
        matched = (record.get('acceptedTaxonKey') or record.get('taxonKey')) == key
    else:
        matched = False
    if not matched or record.get('speciesKey') is None:
        return None, 'requested_taxon_rank_mismatch'
    item, reason = legacy.clean_record(record, record['speciesKey'])
    if item is not None:
        item.update(taxon_key=record.get('taxonKey'),
                    accepted_taxon_key=record.get('acceptedTaxonKey'),
                    species_key=record['speciesKey'],
                    establishment_means=record.get('establishmentMeans'),
                    degree_of_establishment=record.get('degreeOfEstablishment'))
    return item, reason


def acquire_host(taxon: dict, *, fetch: Callable = legacy.get_json) -> dict:
    import requests
    result = {'wcvp_id': taxon['wcvp_id'], 'name': taxon['name'],
              'rank_hint': taxon['rank_hint'], 'status': 'REQUEST_ERROR',
              'quality_records': None, 'points': [], 'request_ledger': [],
              'spatial_filter_applied': False, 'native_range_filter_applied': False,
              'climate_or_genetic_response_used': False}
    rank = taxon['rank_hint']
    if rank not in ALLOWED_RANKS:
        result['status'] = 'TAXONOMIC_RANK_REVIEW'
        return result
    with requests.Session() as session:
        session.headers.update({'User-Agent': 'TTF-host-reference/0.1 (https://github.com/zuizui0223/TTF)'})
        ledger = result['request_ledger']
        try:
            match = fetch(session, 'https://api.gbif.org/v1/species/match',
                {'name': taxon['name'], 'strict': 'true', 'rank': rank, 'kingdom': 'Plantae'}, ledger)
            result['taxonomic_match'] = match
            key = plant_key(match, taxon['name'], rank)
            result['accepted_taxon_key'] = key
            if key is None:
                result['status'] = 'TAXONOMY_NOT_EXACT'
                return result
            query = dict(taxonKey=key, datasetKey=legacy.DATASET, basisOfRecord='HUMAN_OBSERVATION',
                occurrenceStatus='PRESENT', hasCoordinate='true', hasGeospatialIssue='false',
                year='2010,2026', eventDate='2010-01-01,2026-10-18', coordinateUncertaintyInMeters='0,1000')
            count = int(fetch(session, 'https://api.gbif.org/v1/occurrence/search',
                              dict(query, limit=0), ledger)['count'])
            offsets = legacy.window_offsets(count)
            result.update(api_count_initial=count, page_offsets=offsets, api_count_changed=False,
                          retrieval_mode='complete_query' if count <= 900 else 'bounded_3_windows')
            records = []; seen = set(); rejected = Counter(); fetched = 0
            for offset in offsets:
                page = fetch(session, 'https://api.gbif.org/v1/occurrence/search',
                             dict(query, offset=offset, limit=300), ledger)
                if int(page.get('offset', offset)) != offset:
                    raise ValueError('page offset mismatch')
                result['api_count_changed'] |= int(page['count']) != count
                if len(page.get('results', [])) < min(300, max(0, count - offset)):
                    raise ValueError('short page: retain as transport incompleteness')
                for original in page['results']:
                    fetched += 1
                    item, reason = clean_host_record(original, key, rank)
                    if item is None:
                        rejected[reason] += 1; continue
                    if item['gbif_id'] in seen:
                        rejected['duplicate_gbif_id'] += 1; continue
                    seen.add(item['gbif_id']); records.append(item)
            result.update(status='OBSERVATIONS_ACQUIRED', quality_records=len(records), points=records,
                          fetched_records=fetched, rejection_counts=dict(rejected))
        except (requests.RequestException, ValueError, KeyError, RuntimeError, OSError) as exc:
            result.update(status='REQUEST_ERROR', quality_records=None, points=[],
                          error=type(exc).__name__ + ': ' + str(exc)[:400])
    return result


def canonical(value: Any) -> bytes:
    return json.dumps(value, ensure_ascii=False, sort_keys=True, separators=(',', ':'), allow_nan=False).encode()


def validate_result(key: str, result: dict) -> None:
    if result.get('wcvp_id') != key or result.get('status') not in STATUSES:
        raise ValueError('checkpoint identity/status mismatch')
    points = result.get('points')
    if not isinstance(points, list):
        raise ValueError('checkpoint points must be an array')
    if result['status'] == 'OBSERVATIONS_ACQUIRED':
        if result.get('quality_records') != len(points):
            raise ValueError('checkpoint observation count mismatch')
    elif points or result.get('quality_records') is not None:
        raise ValueError('unresolved state cannot be represented as observed zero or retained points')


def save_checkpoint(directory: Path, key: str, binding: dict, result: dict) -> None:
    validate_result(key, result)
    directory = Path(directory); directory.mkdir(parents=True, exist_ok=True)
    payload = {'schema': 'ttf_host_observation_checkpoint_v0.1', 'binding': binding,
               'result': result, 'result_sha256': hashlib.sha256(canonical(result)).hexdigest()}
    fd, temporary = tempfile.mkstemp(prefix='.partial-', dir=directory)
    try:
        with os.fdopen(fd, 'wb') as stream:
            stream.write(canonical(payload) + b'\n'); stream.flush(); os.fsync(stream.fileno())
        os.replace(temporary, directory / (normalize_id(key) + '.json'))
    finally:
        if os.path.exists(temporary): os.unlink(temporary)


def load_checkpoint(directory: Path, key: str, binding: dict) -> dict | None:
    path = Path(directory) / (normalize_id(key) + '.json')
    if not path.exists():
        return None
    payload = json.loads(path.read_text())
    if payload.get('schema') != 'ttf_host_observation_checkpoint_v0.1' or payload.get('binding') != binding:
        raise ValueError('checkpoint source/implementation binding mismatch')
    result = payload['result']
    if hashlib.sha256(canonical(result)).hexdigest() != payload.get('result_sha256'):
        raise ValueError('checkpoint digest mismatch')
    validate_result(key, result)
    return result


def run_queue(plan: list[dict], directory: Path, binding: dict, *, acquire: Callable = acquire_host,
              max_seconds: float = 600, workers: int = 4) -> tuple[list[dict], list[str]]:
    if max_seconds <= 0 or not 1 <= workers <= 4:
        raise ValueError('invalid bounded transport parameters')
    completed = {}; todo = []
    for taxon in plan:
        key = taxon['wcvp_id']
        if key in completed or any(item['wcvp_id'] == key for item in todo):
            raise ValueError('duplicate queue identity')
        previous = load_checkpoint(directory, key, binding)
        if previous is not None: completed[key] = previous
        else: todo.append(taxon)
    pending = iter(todo); start = time.monotonic()
    def worker(taxon):
        result = acquire(taxon)
        save_checkpoint(directory, taxon['wcvp_id'], binding, result)
        return result
    with ThreadPoolExecutor(max_workers=workers) as pool:
        active = {}
        def submit_next():
            if time.monotonic() - start >= max_seconds: return
            taxon = next(pending, None)
            if taxon is not None: active[pool.submit(worker, taxon)] = taxon['wcvp_id']
        for _ in range(workers): submit_next()
        while active:
            finished, _ = wait(active, return_when=FIRST_COMPLETED)
            for future in finished:
                key = active.pop(future); completed[key] = future.result()
            print(json.dumps({'host_taxa_checkpointed': len(completed), 'host_taxa_total': len(plan)}), flush=True)
            for _ in finished: submit_next()
    return ([completed[t['wcvp_id']] for t in plan if t['wcvp_id'] in completed],
            [t['wcvp_id'] for t in plan if t['wcvp_id'] not in completed])


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--manifest', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--max-seconds', type=float, default=600)
    parser.add_argument('--execute', action='store_true', help='Without this, create plan only; no HTTP requests.')
    args = parser.parse_args()
    if hashlib.sha256(Path(legacy.__file__).read_bytes()).hexdigest() != LEGACY_SCRIPT_SHA256:
        raise ValueError('original external-insect acquisition helper SHA mismatch')
    raw = args.manifest.read_bytes()
    if hashlib.sha256(raw).hexdigest() != MANIFEST_SHA256:
        raise ValueError('frozen 1037-association host manifest SHA mismatch')
    with args.manifest.open(newline='', encoding='utf-8') as stream: rows = list(csv.DictReader(stream))
    plan = host_plan(rows)
    if len(rows) != 1037 or len(plan) != 518 or len({r['insect_species'] for r in rows}) != 641:
        raise ValueError('full-panel manifest dimensions changed')
    args.output.mkdir(parents=True, exist_ok=True)
    binding = {'manifest_sha256': MANIFEST_SHA256,
        'host_script_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        'legacy_script_sha256': hashlib.sha256(Path(legacy.__file__).read_bytes()).hexdigest(),
        'plan_sha256': hashlib.sha256(canonical(plan)).hexdigest(), 'dataset': legacy.DATASET,
        'event_date': ['2010-01-01', '2026-10-18']}
    legacy.dump(args.output / 'host_acquisition_plan.json', {'binding': binding, 'taxa': plan})
    if not args.execute:
        print('PLAN_ONLY_518_HOST_TAXA_NO_HTTP_REQUESTS'); return 0
    completed, pending = run_queue(plan, args.output / 'checkpoints', binding, max_seconds=args.max_seconds)
    by_id = {r['wcvp_id']: r for r in completed}
    overview = [{**{k:v for k,v in t.items() if k != 'insect_species'},
                 'status': by_id.get(t['wcvp_id'], {}).get('status', 'NOT_ATTEMPTED'),
                 'quality_records': by_id.get(t['wcvp_id'], {}).get('quality_records')}
                for t in plan]
    with (args.output / 'all_518_host_taxa.csv').open('w', newline='', encoding='utf-8') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(overview[0])); writer.writeheader(); writer.writerows(overview)
    legacy.dump(args.output / 'summary.json', {'status': 'PARTIAL' if pending else 'ALL_TAXA_ATTEMPTED',
        'completed_at': datetime.now(timezone.utc).isoformat(), 'binding': binding,
        'attempted_taxa': len(completed), 'pending_taxa': len(pending), 'pending_ids': pending,
        'status_counts': dict(Counter(x['status'] for x in overview)),
        'quality_records': sum(len(r['points']) for r in completed),
        'spatial_filter_applied': False, 'native_range_filter_applied': False,
        'climate_or_genetic_response_used': False,
        'boundary': 'Modern host observations only. No native or historical occupancy inference; '
                    'each insect-host association still requires its own spatial exclusion and thinning.'})
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
