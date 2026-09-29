#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import json
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

from scripts.acquire_palearctic_lgm_occurrences import (
    in_palearctic,
    load_palearctic_geometry,
)
from ttf.gbif_transport import get_json
from ttf.relational_environment import (
    deterministic_page_offsets,
    filter_and_thin_occurrences,
)
from ttf.relational_external import accepted_gbif_species_match, shard_items


PAGE_WORKERS = 4


def _fetch_page(usage_key: int, offset: int) -> tuple[int, dict]:
    payload = get_json(
        "occurrence/search",
        {
            "taxonKey": usage_key,
            "hasCoordinate": "true",
            "hasGeospatialIssue": "false",
            "occurrenceStatus": "PRESENT",
            "year": "2010,2026",
            "limit": 300,
            "offset": int(offset),
        },
    )
    return int(offset), payload


def fetch_species_palearctic_parallel(
    species: str,
    prepared_realm,
    *,
    minimum_retained: int = 30,
    maximum_retained: int = 200,
) -> tuple[list[dict[str, object]], dict[str, object]]:
    match = get_json("species/match", {"name": species, "strict": "true"})
    accepted, match_rule = accepted_gbif_species_match(species, match)
    usage_key = int(match.get("usageKey") or 0)
    if not accepted or usage_key <= 0:
        return [], {
            "species": species,
            "status": "REJECTED_GBIF_TAXON_MATCH",
            "gbif_usage_key": usage_key,
            "gbif_match_type": match.get("matchType"),
            "gbif_canonical_name": match.get("canonicalName"),
            "gbif_rank": match.get("rank"),
            "match_rule": match_rule,
        }

    count_payload = get_json(
        "occurrence/search",
        {
            "taxonKey": usage_key,
            "hasCoordinate": "true",
            "hasGeospatialIssue": "false",
            "occurrenceStatus": "PRESENT",
            "year": "2010,2026",
            "limit": 1,
        },
    )
    total = int(count_payload.get("count") or 0)
    offsets = deterministic_page_offsets(total, page_size=300, maximum_pages=20)

    # Transport-only acceleration: submit the exact frozen offset set in
    # parallel, then consume results in the original deterministic offset order.
    if offsets:
        with ThreadPoolExecutor(max_workers=min(PAGE_WORKERS, len(offsets))) as pool:
            page_results = list(pool.map(lambda off: _fetch_page(usage_key, off), offsets))
    else:
        page_results = []

    global_records: list[dict[str, object]] = []
    palearctic_records: list[dict[str, object]] = []
    for offset, payload in page_results:
        if offset not in offsets:
            raise RuntimeError("parallel transport returned an unexpected page offset")
        for item in payload.get("results") or []:
            lat = item.get("decimalLatitude")
            lon = item.get("decimalLongitude")
            key = int(item.get("key") or 0)
            if lat is None or lon is None or key <= 0:
                continue
            try:
                lat = float(lat)
                lon = float(lon)
            except (TypeError, ValueError):
                continue
            row = {"source_key": key, "latitude": lat, "longitude": lon}
            global_records.append(row)
            if in_palearctic(prepared_realm, lat, lon):
                palearctic_records.append(row)

    retained = filter_and_thin_occurrences(
        species,
        palearctic_records,
        minimum_distance_km=10.0,
        maximum_retained=maximum_retained,
    )
    status = (
        "PASS_OCCURRENCE_GEOMETRY"
        if len(retained) >= minimum_retained
        else "FAIL_OCCURRENCE_GEOMETRY"
    )
    return retained, {
        "species": species,
        "status": status,
        "gbif_usage_key": usage_key,
        "gbif_total_coordinate_count_2010_2026": total,
        "queried_page_offsets": list(offsets),
        "queried_rows_with_coordinates": len(global_records),
        "queried_rows_in_palearctic": len(palearctic_records),
        "retained_palearctic_after_exact_dedup_thinning_cap": len(retained),
        "gbif_match_type": match.get("matchType"),
        "gbif_canonical_name": match.get("canonicalName"),
        "gbif_rank": match.get("rank"),
        "match_rule": match_rule,
        "transport": {
            "mode": "parallel_frozen_offsets",
            "page_workers": PAGE_WORKERS,
            "scientific_page_offsets_changed": False,
        },
    }


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--census", type=Path, required=True)
    ap.add_argument("--realm-geojson", type=Path, required=True)
    ap.add_argument("--output-csv", type=Path, required=True)
    ap.add_argument("--output-ledger", type=Path, required=True)
    ap.add_argument("--shard-index", type=int, required=True)
    ap.add_argument("--shards", type=int, required=True)
    args = ap.parse_args()

    census = json.loads(args.census.read_text())
    if census.get("schema") != "ttf_genetic_palearctic_lgm_nonlepidoptera_census_v0.2":
        raise RuntimeError("unexpected non-Lepidoptera Palearctic-LGM census")
    if census.get("status") != "PASS_TO_FRESH_V02_GBIF_CENSUS":
        raise RuntimeError("v0.2 census did not authorize external occurrence acquisition")
    if any(bool(v) for v in census["response_firewall"].values()):
        raise RuntimeError("v0.2 genetic response firewall is open")

    names = [str(row["species"]) for row in census["retained_species"]]
    if len(names) != 23 or len(set(names)) != 23:
        raise RuntimeError("expected exact frozen 23-species non-Lepidoptera panel")

    shard = list(shard_items(names, shard_index=args.shard_index, shards=args.shards))
    realm = load_palearctic_geometry(args.realm_geojson)

    output: list[dict[str, object]] = []
    ledger: list[dict[str, object]] = []
    for species in shard:
        try:
            rows, meta = fetch_species_palearctic_parallel(species, realm)
        except Exception as exc:
            rows = []
            meta = {
                "species": species,
                "status": "REQUEST_ERROR",
                "error": f"{type(exc).__name__}: {exc}",
                "transport": {
                    "mode": "parallel_frozen_offsets",
                    "page_workers": PAGE_WORKERS,
                    "scientific_page_offsets_changed": False,
                },
            }
        for row in rows:
            row["role"] = "both"
        output.extend(rows)
        meta["role"] = "both"
        ledger.append(meta)
        print(json.dumps({
            "species": species,
            "status": meta["status"],
            "retained": meta.get("retained_palearctic_after_exact_dedup_thinning_cap", 0),
        }, sort_keys=True), flush=True)

    args.output_csv.parent.mkdir(parents=True, exist_ok=True)
    fields = [
        "species", "role", "source_key", "latitude", "longitude",
        "priority_rank", "priority_sha256",
    ]
    with args.output_csv.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=fields)
        writer.writeheader()
        writer.writerows(output)

    payload = {
        "schema": "ttf_genetic_palearctic_lgm_occurrence_shard_v0.2",
        "shard_index": int(args.shard_index),
        "shards": int(args.shards),
        "species_requested": len(shard),
        "retained_occurrence_rows": len(output),
        "species": ledger,
        "transport": {
            "mode": "parallel_frozen_offsets",
            "page_workers": PAGE_WORKERS,
            "scientific_page_offsets_changed": False,
        },
        "response_firewall": {
            "species_level_genetic_scores_used": False,
            "pairwise_v02_T_st_computed": False,
            "beta_LGM_computed": False,
        },
    }
    args.output_ledger.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
