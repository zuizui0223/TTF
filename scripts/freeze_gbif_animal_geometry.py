#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import time
from pathlib import Path
from urllib.parse import urlencode
from urllib.request import Request, urlopen

EARTH_RADIUS_KM = 6371.0088
GBIF = "https://api.gbif.org/v1"
USER_AGENT = "ttf-gate-i-geometry/0.1 (https://github.com/zuizui0223/TTF)"


def get_json(path: str, params: dict[str, object], retries: int = 4) -> dict:
    url = f"{GBIF}/{path}?{urlencode(params)}"
    error: Exception | None = None
    for attempt in range(retries):
        try:
            request = Request(
                url,
                headers={"User-Agent": USER_AGENT, "Accept": "application/json"},
            )
            with urlopen(request, timeout=60) as response:
                return json.loads(response.read().decode("utf-8"))
        except Exception as exc:  # network boundary
            error = exc
            time.sleep(2**attempt)
    raise RuntimeError(f"GBIF request failed: {url}") from error


def polygon(row: dict[str, str]) -> str:
    x0, y0 = float(row["min_lon"]), float(row["min_lat"])
    x1, y1 = float(row["max_lon"]), float(row["max_lat"])
    return f"POLYGON(({x0} {y0},{x1} {y0},{x1} {y1},{x0} {y1},{x0} {y0}))"


def hash_priority(seed: int, species: str, occurrence_key: int) -> str:
    return hashlib.sha256(
        f"{int(seed)}|{species}|{int(occurrence_key)}".encode("utf-8")
    ).hexdigest()


def xyz_km(latitude: float, longitude: float) -> tuple[float, float, float]:
    lat = math.radians(latitude)
    lon = math.radians(longitude)
    cos_lat = math.cos(lat)
    return (
        EARTH_RADIUS_KM * cos_lat * math.cos(lon),
        EARTH_RADIUS_KM * cos_lat * math.sin(lon),
        EARTH_RADIUS_KM * math.sin(lat),
    )


def normalized_name(value: object) -> str:
    return " ".join(str(value or "").strip().split()).casefold()


def accepted_taxon_match(
    species: str,
    match: dict[str, object],
    *,
    allow_canonical_fuzzy: bool,
) -> tuple[bool, str]:
    """Return whether a GBIF match is admissible and the frozen acceptance rule.

    The historical Gate-I-A default remains unchanged: only EXACT/CONFIDENCE
    matches are accepted. A fresh prospective panel may explicitly opt into a
    narrow fallback for GBIF ``FUZZY`` records only when GBIF's canonical name is
    exactly the requested species name and the matched rank is SPECIES. This
    handles matcher metadata quirks without silently substituting a different
    taxon.
    """
    match_type = str(match.get("matchType") or "").upper()
    if match_type in {"EXACT", "CONFIDENCE"}:
        return True, "exact_or_confidence"
    if not allow_canonical_fuzzy or match_type != "FUZZY":
        return False, "rejected_match_type"

    canonical = normalized_name(match.get("canonicalName"))
    requested = normalized_name(species)
    rank = str(match.get("rank") or "").upper()
    if canonical == requested and rank == "SPECIES":
        return True, "canonical_name_exact_species_rank_fuzzy"
    return False, "rejected_fuzzy_without_exact_canonical_species_match"


def fetch_species(
    row: dict[str, str],
    *,
    fetch_limit: int,
    cap: int,
    seed: int,
    allow_canonical_fuzzy: bool = False,
) -> tuple[list[dict[str, object]], dict[str, object]]:
    species = row["scientific_name"].strip()
    match = get_json("species/match", {"name": species, "strict": "true"})
    usage_key = int(match.get("usageKey") or 0)
    if not usage_key:
        raise RuntimeError(f"GBIF exact taxon match failed for {species}")
    match_type = str(match.get("matchType") or "")
    accepted, acceptance_rule = accepted_taxon_match(
        species,
        match,
        allow_canonical_fuzzy=allow_canonical_fuzzy,
    )
    if not accepted:
        raise RuntimeError(
            "GBIF taxon match for "
            f"{species} is not admissible: matchType={match_type!r}, "
            f"canonicalName={match.get('canonicalName')!r}, rank={match.get('rank')!r}, "
            f"rule={acceptance_rule!r}"
        )

    records: list[dict[str, object]] = []
    offset = 0
    while len(records) < fetch_limit:
        payload = get_json(
            "occurrence/search",
            {
                "taxon_key": usage_key,
                "geometry": polygon(row),
                "has_coordinate": "true",
                "has_geospatial_issue": "false",
                "occurrence_status": "present",
                "limit": min(300, fetch_limit - len(records)),
                "offset": offset,
            },
        )
        batch = payload.get("results", [])
        if not batch:
            break
        for item in batch:
            latitude = item.get("decimalLatitude")
            longitude = item.get("decimalLongitude")
            key = int(item.get("key") or 0)
            if latitude is None or longitude is None or key <= 0:
                continue
            lat = float(latitude)
            lon = float(longitude)
            if not (-90 <= lat <= 90 and -180 <= lon <= 180):
                continue
            records.append(
                {
                    "source_key": key,
                    "latitude": lat,
                    "longitude": lon,
                }
            )
        offset += len(batch)
        if payload.get("endOfRecords", False):
            break

    # The benchmark keeps one record per exact coordinate. Dense nearby records
    # remain; only literal duplicates are removed. A hash-priority cap makes the
    # retained frame deterministic given the fetched GBIF keys and frozen seed.
    by_coordinate: dict[tuple[float, float], dict[str, object]] = {}
    for record in records:
        coordinate = (float(record["latitude"]), float(record["longitude"]))
        current = by_coordinate.get(coordinate)
        if current is None or int(record["source_key"]) < int(current["source_key"]):
            by_coordinate[coordinate] = record
    unique = list(by_coordinate.values())
    unique.sort(
        key=lambda record: (
            hash_priority(seed, species, int(record["source_key"])),
            int(record["source_key"]),
        )
    )
    retained = unique[:cap]
    for record in retained:
        x, y, z = xyz_km(float(record["latitude"]), float(record["longitude"]))
        record.update(species=species, x_km=x, y_km=y, z_km=z)

    keys = sorted(int(record["source_key"]) for record in retained)
    ledger = {
        "species": species,
        "group": row.get("group", ""),
        "gbif_usage_key": usage_key,
        "gbif_match_type": match_type,
        "gbif_matched_canonical_name": match.get("canonicalName"),
        "gbif_matched_rank": match.get("rank"),
        "gbif_match_acceptance_rule": acceptance_rule,
        "query_geometry": polygon(row),
        "fetched_coordinate_rows": len(records),
        "unique_exact_coordinates": len(unique),
        "retained_records": len(retained),
        "retained_source_keys_sha256": hashlib.sha256(
            json.dumps(keys, separators=(",", ":")).encode("utf-8")
        ).hexdigest(),
    }
    return retained, ledger


def canonical_sha256(value: object) -> str:
    return hashlib.sha256(
        json.dumps(
            value,
            sort_keys=True,
            separators=(",", ":"),
            ensure_ascii=True,
        ).encode("utf-8")
    ).hexdigest()


def main() -> int:
    parser = argparse.ArgumentParser(
        description=(
            "Freeze a non-flower GBIF animal occurrence sampling geometry for "
            "TTF Gate-I-A external stress qualification."
        )
    )
    parser.add_argument(
        "--manifest",
        type=Path,
        default=Path("benchmarks/gate_i_animal_geometry_manifest.csv"),
    )
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--ledger", type=Path, required=True)
    parser.add_argument("--fetch-limit", type=int, default=1200)
    parser.add_argument("--cap-per-species", type=int, default=120)
    parser.add_argument("--minimum-records", type=int, default=80)
    parser.add_argument("--seed", type=int, default=20260907)
    parser.add_argument(
        "--allow-canonical-fuzzy",
        action="store_true",
        help=(
            "Also accept GBIF FUZZY matches only when canonicalName exactly equals "
            "the requested name and rank is SPECIES. Default false preserves the "
            "historical Gate-I-A acquisition rule."
        ),
    )
    args = parser.parse_args()
    if args.fetch_limit < args.cap_per_species:
        raise ValueError("fetch-limit must be >= cap-per-species")
    if not 8 <= args.minimum_records <= args.cap_per_species:
        raise ValueError("minimum-records must lie in [8, cap-per-species]")

    declarations = list(csv.DictReader(args.manifest.open(encoding="utf-8")))
    if len(declarations) < 12:
        raise ValueError("external geometry benchmark requires at least 12 taxa")
    names = [row["scientific_name"].strip() for row in declarations]
    if any(not name for name in names) or len(set(names)) != len(names):
        raise ValueError("manifest scientific names must be unique and non-empty")

    all_records: list[dict[str, object]] = []
    species_ledgers: list[dict[str, object]] = []
    for row in declarations:
        records, ledger = fetch_species(
            row,
            fetch_limit=args.fetch_limit,
            cap=args.cap_per_species,
            seed=args.seed,
            allow_canonical_fuzzy=args.allow_canonical_fuzzy,
        )
        species_ledgers.append(ledger)
        if len(records) < args.minimum_records:
            raise RuntimeError(
                f"{row['scientific_name']} retained {len(records)} records; "
                f"minimum is {args.minimum_records}"
            )
        all_records.extend(records)

    all_records.sort(
        key=lambda record: (str(record["species"]), int(record["source_key"]))
    )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", newline="", encoding="utf-8") as handle:
        fields = [
            "species",
            "x_km",
            "y_km",
            "z_km",
            "source_key",
            "latitude",
            "longitude",
        ]
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(all_records)

    source_key_pairs = [
        [str(record["species"]), int(record["source_key"])] for record in all_records
    ]
    payload = {
        "schema": "ttf_gate_i_a_gbif_animal_geometry_v0.1",
        "purpose": "external_nonflower_sampling_geometry_stress_only",
        "provider": "GBIF occurrence/search",
        "manifest": str(args.manifest),
        "manifest_sha256": hashlib.sha256(args.manifest.read_bytes()).hexdigest(),
        "seed": args.seed,
        "fetch_limit": args.fetch_limit,
        "cap_per_species": args.cap_per_species,
        "minimum_records": args.minimum_records,
        "allow_canonical_fuzzy": bool(args.allow_canonical_fuzzy),
        "species_count": len(species_ledgers),
        "record_count": len(all_records),
        "source_key_pairs_sha256": canonical_sha256(source_key_pairs),
        "species": species_ledgers,
        "claim_boundary": (
            "This freezes an external animal occurrence geometry stress fixture. "
            "It does not replace Gate-I-B qualification on the intended empirical "
            "TTF sampling geometry."
        ),
    }
    payload["ledger_sha256"] = canonical_sha256(payload)
    args.ledger.parent.mkdir(parents=True, exist_ok=True)
    args.ledger.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    print(json.dumps({
        "species_count": payload["species_count"],
        "record_count": payload["record_count"],
        "ledger_sha256": payload["ledger_sha256"],
    }, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
