#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import json
import time
from pathlib import Path
from urllib.error import HTTPError
from urllib.parse import urlencode
from urllib.request import Request, urlopen

from ttf.relational_environment import (
    deterministic_page_offsets,
    filter_and_thin_occurrences,
)
from ttf.relational_external import accepted_gbif_species_match, shard_items

GBIF = "https://api.gbif.org/v1"
USER_AGENT = "ttf-palearctic-insect-lgm/0.1"


def get_json(path: str, params: dict[str, object], retries: int = 8) -> dict:
    url = f"{GBIF}/{path}?{urlencode(params)}"
    error: Exception | None = None
    for attempt in range(retries):
        try:
            req = Request(url, headers={"User-Agent": USER_AGENT, "Accept": "application/json"})
            with urlopen(req, timeout=90) as response:
                out = json.loads(response.read().decode("utf-8"))
            time.sleep(0.05)
            return out
        except HTTPError as exc:
            error = exc
            if exc.code not in {429, 500, 502, 503, 504}:
                break
            retry_after = exc.headers.get("Retry-After")
            try:
                wait = float(retry_after) if retry_after else float(2**attempt)
            except (TypeError, ValueError):
                wait = float(2**attempt)
            time.sleep(min(120.0, max(1.0, wait)))
        except Exception as exc:
            error = exc
            time.sleep(min(60.0, float(2**attempt)))
    raise RuntimeError(f"GBIF request failed: {url}") from error


def fetch_species(species: str) -> tuple[list[dict[str, object]], dict[str, object]]:
    match = get_json("species/match", {"name": species, "strict": "true"})
    accepted, match_rule = accepted_gbif_species_match(species, match)
    usage_key = int(match.get("usageKey") or 0)
    taxonomy = {
        "kingdom": match.get("kingdom"),
        "phylum": match.get("phylum"),
        "class": match.get("class"),
        "order": match.get("order"),
        "family": match.get("family"),
        "genus": match.get("genus"),
    }
    if not accepted or usage_key <= 0:
        return [], {
            "species": species,
            "status": "REJECTED_GBIF_TAXON_MATCH",
            "gbif_usage_key": usage_key,
            "gbif_match_type": match.get("matchType"),
            "gbif_canonical_name": match.get("canonicalName"),
            "gbif_rank": match.get("rank"),
            "match_rule": match_rule,
            "taxonomy": taxonomy,
        }

    count = get_json(
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
    total = int(count.get("count") or 0)
    offsets = deterministic_page_offsets(total, page_size=300, maximum_pages=20)
    raw: list[dict[str, object]] = []
    for offset in offsets:
        page = get_json(
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
        for item in page.get("results") or []:
            lat = item.get("decimalLatitude")
            lon = item.get("decimalLongitude")
            key = int(item.get("key") or 0)
            if lat is None or lon is None or key <= 0:
                continue
            raw.append({"source_key": key, "latitude": lat, "longitude": lon})

    retained = filter_and_thin_occurrences(
        species,
        raw,
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
        "queried_rows_with_coordinates": len(raw),
        "retained_after_exact_dedup_thinning_cap": len(retained),
        "gbif_match_type": match.get("matchType"),
        "gbif_canonical_name": match.get("canonicalName"),
        "gbif_rank": match.get("rank"),
        "match_rule": match_rule,
        "taxonomy": taxonomy,
    }


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--candidates", type=Path, required=True)
    ap.add_argument("--output-csv", type=Path, required=True)
    ap.add_argument("--output-ledger", type=Path, required=True)
    ap.add_argument("--shard-index", type=int, required=True)
    ap.add_argument("--shards", type=int, required=True)
    args = ap.parse_args()

    candidate = json.loads(args.candidates.read_text())
    if candidate.get("schema") != "ttf_palearctic_insect_lgm_taxonomic_candidates_v0.1":
        raise RuntimeError("unexpected candidate manifest")
    names = list(map(str, candidate["species"]))
    if len(names) != 143 or len(set(names)) != 143:
        raise RuntimeError("taxonomic parent must contain exactly 143 unique species")
    shard = list(shard_items(names, shard_index=args.shard_index, shards=args.shards))

    retained_rows: list[dict[str, object]] = []
    ledgers: list[dict[str, object]] = []
    for species in shard:
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

    args.output_csv.parent.mkdir(parents=True, exist_ok=True)
    fields = ["species", "source_key", "latitude", "longitude", "priority_rank", "priority_sha256"]
    with args.output_csv.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(retained_rows)

    payload = {
        "schema": "ttf_palearctic_insect_lgm_gbif_shard_v0.1",
        "shard_index": args.shard_index,
        "shards": args.shards,
        "species_requested": len(shard),
        "species": ledgers,
        "response_firewall": {
            "species_level_Phase4_scores_used": False,
            "pairwise_genetic_distances_used": False,
            "subpanel_T_st_computed": False,
            "subpanel_beta_LGM_computed": False,
        },
    }
    args.output_ledger.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
