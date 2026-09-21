#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import json
import time
from pathlib import Path
from urllib.parse import urlencode
from urllib.request import Request, urlopen

from ttf.relational_environment import deterministic_page_offsets, filter_and_thin_occurrences
from ttf.relational_external import accepted_gbif_species_match, shard_items

GBIF = "https://api.gbif.org/v1"
USER_AGENT = "ttf-relational-environment/0.2 (https://github.com/zuizui0223/TTF)"


def get_json(path: str, params: dict[str, object], retries: int = 5) -> dict:
    url = f"{GBIF}/{path}?{urlencode(params)}"
    error: Exception | None = None
    for attempt in range(retries):
        try:
            request = Request(url, headers={"User-Agent": USER_AGENT, "Accept": "application/json"})
            with urlopen(request, timeout=60) as response:
                return json.loads(response.read().decode("utf-8"))
        except Exception as exc:
            error = exc
            time.sleep(min(30, 2**attempt))
    raise RuntimeError(f"GBIF request failed: {url}") from error


def fetch_species(species: str) -> tuple[list[dict[str, object]], dict[str, object]]:
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
            "taxon_key": usage_key,
            "has_coordinate": "true",
            "has_geospatial_issue": "false",
            "occurrence_status": "present",
            "year": "2010,2026",
            "limit": 1,
        },
    )
    total = int(count_payload.get("count") or 0)
    offsets = deterministic_page_offsets(total, page_size=300, maximum_pages=20)
    records: list[dict[str, object]] = []
    for offset in offsets:
        payload = get_json(
            "occurrence/search",
            {
                "taxon_key": usage_key,
                "has_coordinate": "true",
                "has_geospatial_issue": "false",
                "occurrence_status": "present",
                "year": "2010,2026",
                "limit": 300,
                "offset": int(offset),
            },
        )
        for item in payload.get("results") or []:
            lat = item.get("decimalLatitude")
            lon = item.get("decimalLongitude")
            key = int(item.get("key") or 0)
            if lat is None or lon is None or key <= 0:
                continue
            records.append({"source_key": key, "latitude": lat, "longitude": lon})

    retained = filter_and_thin_occurrences(
        species,
        records,
        minimum_distance_km=10.0,
        maximum_retained=200,
    )
    status = "PASS_OCCURRENCE_GEOMETRY" if len(retained) >= 30 else "FAIL_OCCURRENCE_GEOMETRY"
    return retained, {
        "species": species,
        "status": status,
        "gbif_usage_key": usage_key,
        "gbif_total_coordinate_count_2010_2026": total,
        "queried_page_offsets": list(offsets),
        "queried_rows_with_coordinates": len(records),
        "retained_after_exact_dedup_thinning_cap": len(retained),
        "gbif_match_type": match.get("matchType"),
        "gbif_canonical_name": match.get("canonicalName"),
        "gbif_rank": match.get("rank"),
        "match_rule": match_rule,
    }


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--candidates", type=Path, required=True)
    ap.add_argument("--output-csv", type=Path, required=True)
    ap.add_argument("--output-ledger", type=Path, required=True)
    ap.add_argument("--shard-index", type=int, default=0)
    ap.add_argument("--shards", type=int, default=1)
    args = ap.parse_args()

    rows = list(csv.DictReader(args.candidates.open(encoding="utf-8")))
    names = [str(row["species"]).strip() for row in rows]
    if len(names) != 1000 or len(set(names)) != 1000:
        raise RuntimeError("expected exact frozen 1000-species candidate set")
    shard = list(shard_items(names, shard_index=args.shard_index, shards=args.shards))

    retained_rows: list[dict[str, object]] = []
    ledgers: list[dict[str, object]] = []
    for i, species in enumerate(shard, start=1):
        try:
            retained, ledger = fetch_species(species)
        except Exception as exc:
            retained, ledger = [], {
                "species": species,
                "status": "REQUEST_ERROR",
                "error": f"{type(exc).__name__}: {exc}",
            }
        retained_rows.extend(retained)
        ledgers.append(ledger)
        if i % 5 == 0:
            print(json.dumps({"shard": args.shard_index, "completed": i, "total": len(shard)}))

    args.output_csv.parent.mkdir(parents=True, exist_ok=True)
    fields = ["species", "source_key", "latitude", "longitude", "priority_rank", "priority_sha256"]
    with args.output_csv.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(retained_rows)

    payload = {
        "schema": "ttf_relational_environment_occurrence_shard_v0.2",
        "shard_index": args.shard_index,
        "shards": args.shards,
        "species_requested": len(shard),
        "retained_occurrence_rows": len(retained_rows),
        "species": ledgers,
        "response_firewall": {
            "study_B_sequence_identity_opened": False,
            "study_B_pairwise_genetic_distances_opened": False,
            "study_B_T_st_computed": False,
            "study_B_beta_R_computed": False,
        },
    }
    args.output_ledger.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
