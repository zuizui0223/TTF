#!/usr/bin/env python3
from __future__ import annotations

import argparse
import hashlib
import json
import time
from pathlib import Path
from urllib.parse import urlencode
from urllib.request import Request, urlopen

from ttf.relational_external import (
    accepted_gbif_species_match,
    deterministic_greedy_thin,
    shard_items,
)

GBIF = "https://api.gbif.org/v1"
USER_AGENT = "ttf-relational-feasibility/0.1 (https://github.com/zuizui0223/TTF)"


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


def fetch_species(
    species: str,
    *,
    year_range: str,
    thin_km: float,
    required_thinned: int,
    fetch_cap: int,
) -> dict:
    match = get_json("species/match", {"name": species, "strict": "true"})
    accepted, rule = accepted_gbif_species_match(species, match)
    if not accepted:
        return {
            "species": species,
            "status": "taxon_match_rejected",
            "match_type": match.get("matchType"),
            "canonical_name": match.get("canonicalName"),
            "match_rank": match.get("rank"),
            "match_rule": rule,
            "fetched": 0,
            "unique_exact_coordinates": 0,
            "thinned_records": 0,
        }

    usage_key = int(match["usageKey"])
    records: list[dict[str, object]] = []
    offset = 0
    end_of_records = False
    total_count = None
    while len(records) < fetch_cap and not end_of_records:
        limit = min(300, fetch_cap - len(records))
        payload = get_json(
            "occurrence/search",
            {
                "taxon_key": usage_key,
                "has_coordinate": "true",
                "has_geospatial_issue": "false",
                "occurrence_status": "present",
                "year": year_range,
                "limit": limit,
                "offset": offset,
            },
        )
        if total_count is None:
            total_count = int(payload.get("count") or 0)
        batch = payload.get("results") or []
        if not batch:
            end_of_records = True
            break
        for item in batch:
            lat = item.get("decimalLatitude")
            lon = item.get("decimalLongitude")
            key = int(item.get("key") or 0)
            if lat is None or lon is None or key <= 0:
                continue
            records.append({"source_key": key, "latitude": lat, "longitude": lon})
        offset += len(batch)
        end_of_records = bool(payload.get("endOfRecords", False))
        # Once a strong lower bound is reached, one additional page is not
        # required: this census only asks whether >= required_thinned is feasible.
        thinned = deterministic_greedy_thin(records, minimum_distance_km=thin_km)
        if len(thinned) >= required_thinned:
            break

    thinned = deterministic_greedy_thin(records, minimum_distance_km=thin_km)
    unique_coords = {
        (float(r["latitude"]), float(r["longitude"]))
        for r in records
        if r.get("latitude") is not None and r.get("longitude") is not None
    }
    if len(thinned) >= required_thinned:
        status = "usable_lower_bound"
    elif end_of_records:
        status = "insufficient_complete_search"
    else:
        status = "censored_below_threshold_at_fetch_cap"
    return {
        "species": species,
        "status": status,
        "gbif_usage_key": usage_key,
        "match_type": match.get("matchType"),
        "canonical_name": match.get("canonicalName"),
        "match_rank": match.get("rank"),
        "match_rule": rule,
        "gbif_total_coordinate_count": total_count,
        "fetched": len(records),
        "unique_exact_coordinates": len(unique_coords),
        "thinned_records": len(thinned),
        "retained_source_keys_sha256": hashlib.sha256(
            json.dumps(sorted(int(r["source_key"]) for r in thinned), separators=(",", ":")).encode()
        ).hexdigest(),
    }


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--census", type=Path, required=True)
    ap.add_argument("--output", type=Path, required=True)
    ap.add_argument("--shard-index", type=int, default=0)
    ap.add_argument("--shards", type=int, default=1)
    ap.add_argument("--year-range", default="2010,2026")
    ap.add_argument("--thin-km", type=float, default=10.0)
    ap.add_argument("--required-thinned", type=int, default=30)
    ap.add_argument("--fetch-cap", type=int, default=1200)
    args = ap.parse_args()

    source = json.loads(args.census.read_text())
    if source.get("schema") != "ttf_relational_fresh_species_census_v0.1":
        raise RuntimeError("unexpected fresh-species census schema")
    firewall = source.get("outcome_firewall") or {}
    if any(bool(v) for v in firewall.values()):
        raise RuntimeError("fresh-species census outcome firewall is open")

    selected = list(source.get("selected_species") or [])
    names = [str(row["species"]) for row in selected]
    shard = shard_items(names, shard_index=args.shard_index, shards=args.shards)
    results = []
    for i, species in enumerate(shard, start=1):
        try:
            results.append(
                fetch_species(
                    species,
                    year_range=args.year_range,
                    thin_km=args.thin_km,
                    required_thinned=args.required_thinned,
                    fetch_cap=args.fetch_cap,
                )
            )
        except Exception as exc:
            results.append({"species": species, "status": "request_error", "error": str(exc)})
        if i % 10 == 0:
            print(json.dumps({"shard": args.shard_index, "completed": i, "total": len(shard)}))

    payload = {
        "schema": "ttf_relational_gbif_occurrence_feasibility_shard_v0.1",
        "source_census_digest_sha256": source["selected_species_digest_sha256"],
        "shard_index": args.shard_index,
        "shards": args.shards,
        "year_range": args.year_range,
        "thin_km": args.thin_km,
        "required_thinned": args.required_thinned,
        "fetch_cap": args.fetch_cap,
        "results": results,
        "firewall": {
            "nucleotide_identity_opened": False,
            "pairwise_genetic_distances_opened": False,
            "T_st_computed": False,
            "ecology_genetics_association_computed": False,
        },
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
