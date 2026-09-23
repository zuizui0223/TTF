from __future__ import annotations

import math
from typing import Iterable, Mapping


EARTH_RADIUS_KM = 6371.0088


def normalized_name(value: object) -> str:
    return " ".join(str(value or "").strip().split()).casefold()


def accepted_gbif_species_match(
    requested_species: str,
    match: Mapping[str, object],
) -> tuple[bool, str]:
    """Strict response-blind GBIF taxon acceptance for relational feasibility."""
    usage_key = int(match.get("usageKey") or 0)
    if usage_key <= 0:
        return False, "missing_usage_key"
    match_type = str(match.get("matchType") or "").upper()
    rank = str(match.get("rank") or "").upper()
    canonical = normalized_name(match.get("canonicalName"))
    requested = normalized_name(requested_species)
    if match_type in {"EXACT", "CONFIDENCE"} and rank == "SPECIES" and canonical == requested:
        return True, "exact_canonical_species"
    return False, "rejected_nonexact_species_match"


def haversine_km(a: tuple[float, float], b: tuple[float, float]) -> float:
    lat1, lon1 = map(math.radians, a)
    lat2, lon2 = map(math.radians, b)
    dlat = lat2 - lat1
    dlon = lon2 - lon1
    h = math.sin(dlat / 2.0) ** 2 + math.cos(lat1) * math.cos(lat2) * math.sin(dlon / 2.0) ** 2
    return 2.0 * EARTH_RADIUS_KM * math.asin(min(1.0, math.sqrt(h)))


def deterministic_greedy_thin(
    records: Iterable[Mapping[str, object]],
    *,
    minimum_distance_km: float,
) -> list[dict[str, object]]:
    """Greedy spatial thinning after deterministic GBIF-key ordering.

    This is only a feasibility lower-bound screen. It is not the final niche
    sampling representation used by the confirmatory relational analysis.
    """
    if minimum_distance_km <= 0:
        raise ValueError("minimum_distance_km must be positive")
    cleaned: dict[tuple[float, float], dict[str, object]] = {}
    for row in records:
        try:
            key = int(row["source_key"])
            lat = float(row["latitude"])
            lon = float(row["longitude"])
        except (KeyError, TypeError, ValueError):
            continue
        if key <= 0 or not (-90.0 <= lat <= 90.0 and -180.0 <= lon <= 180.0):
            continue
        coord = (lat, lon)
        current = cleaned.get(coord)
        if current is None or key < int(current["source_key"]):
            cleaned[coord] = {"source_key": key, "latitude": lat, "longitude": lon}
    ordered = sorted(cleaned.values(), key=lambda r: int(r["source_key"]))
    retained: list[dict[str, object]] = []
    for row in ordered:
        coord = (float(row["latitude"]), float(row["longitude"]))
        if all(
            haversine_km(coord, (float(prev["latitude"]), float(prev["longitude"])))
            >= minimum_distance_km
            for prev in retained
        ):
            retained.append(row)
    return retained


def shard_items(items: list[object], *, shard_index: int, shards: int) -> list[object]:
    if shards < 1 or not 0 <= shard_index < shards:
        raise ValueError("invalid shard specification")
    return [item for i, item in enumerate(items) if i % shards == shard_index]
