from __future__ import annotations

from dataclasses import dataclass
import hashlib
import json
import math
from pathlib import Path
from typing import Iterable, Mapping, Sequence


@dataclass(frozen=True)
class SpeciesRealmSummary:
    species: str
    class_name: str
    order: str
    family: str
    total_localities: int
    palearctic_localities: int
    palearctic_fraction: float
    terrestrial_core: bool
    realm_eligible: bool
    biologically_eligible: bool


def _same_point(a: Sequence[float], b: Sequence[float], tol: float = 1e-12) -> bool:
    return abs(float(a[0]) - float(b[0])) <= tol and abs(float(a[1]) - float(b[1])) <= tol


def _point_on_segment(
    x: float,
    y: float,
    a: Sequence[float],
    b: Sequence[float],
    *,
    tol: float = 1e-10,
) -> bool:
    ax, ay = float(a[0]), float(a[1])
    bx, by = float(b[0]), float(b[1])
    cross = (x - ax) * (by - ay) - (y - ay) * (bx - ax)
    scale = max(1.0, abs(bx - ax), abs(by - ay))
    if abs(cross) > tol * scale:
        return False
    dot = (x - ax) * (x - bx) + (y - ay) * (y - by)
    return dot <= tol * scale * scale


def _ring_contains_or_boundary(
    lon: float,
    lat: float,
    ring: Sequence[Sequence[float]],
) -> tuple[bool, bool]:
    """Return (inside, on_boundary) for one GeoJSON linear ring."""
    if len(ring) < 4:
        raise ValueError("GeoJSON ring must contain at least four coordinates")
    points = list(ring)
    if not _same_point(points[0], points[-1]):
        points.append(points[0])

    inside = False
    for a, b in zip(points[:-1], points[1:]):
        if _point_on_segment(lon, lat, a, b):
            return True, True
        x1, y1 = float(a[0]), float(a[1])
        x2, y2 = float(b[0]), float(b[1])
        if (y1 > lat) == (y2 > lat):
            continue
        x_intersect = x1 + (lat - y1) * (x2 - x1) / (y2 - y1)
        if lon < x_intersect:
            inside = not inside
    return inside, False


def _polygon_contains_or_boundary(
    lon: float,
    lat: float,
    polygon: Sequence[Sequence[Sequence[float]]],
) -> bool:
    if not polygon:
        return False
    outer_inside, outer_boundary = _ring_contains_or_boundary(lon, lat, polygon[0])
    if outer_boundary:
        return True
    if not outer_inside:
        return False
    for hole in polygon[1:]:
        hole_inside, hole_boundary = _ring_contains_or_boundary(lon, lat, hole)
        if hole_boundary:
            return True
        if hole_inside:
            return False
    return True


def geometry_contains_or_boundary(
    geometry: Mapping[str, object],
    *,
    lon: float,
    lat: float,
) -> bool:
    kind = str(geometry.get("type", ""))
    coordinates = geometry.get("coordinates")
    if kind == "Polygon":
        if not isinstance(coordinates, list):
            raise ValueError("Polygon coordinates must be a list")
        return _polygon_contains_or_boundary(lon, lat, coordinates)
    if kind == "MultiPolygon":
        if not isinstance(coordinates, list):
            raise ValueError("MultiPolygon coordinates must be a list")
        return any(_polygon_contains_or_boundary(lon, lat, polygon) for polygon in coordinates)
    raise ValueError(f"unsupported GeoJSON geometry type: {kind!r}")


def select_realm_geometries(
    geojson: Mapping[str, object],
    *,
    realm: str = "Palearctic",
    realmcode: str = "PA",
) -> tuple[Mapping[str, object], ...]:
    if geojson.get("type") != "FeatureCollection":
        raise ValueError("authoritative realm input must be a GeoJSON FeatureCollection")
    features = geojson.get("features")
    if not isinstance(features, list):
        raise ValueError("GeoJSON FeatureCollection lacks features")

    selected: list[Mapping[str, object]] = []
    for feature in features:
        if not isinstance(feature, Mapping):
            continue
        props = feature.get("properties")
        geom = feature.get("geometry")
        if not isinstance(props, Mapping) or not isinstance(geom, Mapping):
            continue
        code = str(props.get("realmcode", props.get("REALMCODE", ""))).strip()
        name = str(props.get("realm", props.get("REALM", ""))).strip()
        if code == realmcode or name.casefold() == realm.casefold():
            selected.append(geom)
    if not selected:
        raise ValueError(f"no GeoJSON feature matched realm={realm!r} realmcode={realmcode!r}")
    return tuple(selected)


def point_in_any_geometry(
    lon: float,
    lat: float,
    geometries: Sequence[Mapping[str, object]],
) -> bool:
    if not math.isfinite(float(lon)) or not math.isfinite(float(lat)):
        raise ValueError("locality coordinates must be finite")
    return any(
        geometry_contains_or_boundary(geometry, lon=float(lon), lat=float(lat))
        for geometry in geometries
    )


def terrestrial_core(
    *,
    class_name: str,
    order: str,
    family: str,
    allowed_orders: Iterable[str],
    excluded_coleoptera_families: Iterable[str],
) -> bool:
    if str(class_name) != "Insecta":
        return False
    if str(order) not in set(map(str, allowed_orders)):
        return False
    if str(order) == "Coleoptera" and str(family) in set(map(str, excluded_coleoptera_families)):
        return False
    return True


def summarize_species_realm(
    rows: Iterable[Mapping[str, object]],
    *,
    geometries: Sequence[Mapping[str, object]],
    allowed_orders: Iterable[str],
    excluded_coleoptera_families: Iterable[str],
    minimum_fraction: float,
    minimum_localities: int,
) -> tuple[SpeciesRealmSummary, ...]:
    if not (0.0 <= float(minimum_fraction) <= 1.0):
        raise ValueError("minimum_fraction must lie in [0, 1]")
    if int(minimum_localities) < 1:
        raise ValueError("minimum_localities must be positive")

    grouped: dict[str, list[Mapping[str, object]]] = {}
    for row in rows:
        species = str(row["species"]).strip()
        if not species:
            raise ValueError("blank species")
        grouped.setdefault(species, []).append(row)

    summaries: list[SpeciesRealmSummary] = []
    for species in sorted(grouped):
        records = grouped[species]
        taxonomy = {
            (str(row["class"]).strip(), str(row["order"]).strip(), str(row["family"]).strip())
            for row in records
        }
        if len(taxonomy) != 1:
            raise ValueError(f"taxonomy varies within species: {species}")
        class_name, order, family = next(iter(taxonomy))

        points = {
            (float(row["longitude"]), float(row["latitude"]))
            for row in records
        }
        if len(points) != len(records):
            raise ValueError(f"duplicate response-blind localities found for species: {species}")
        pal = sum(
            point_in_any_geometry(lon, lat, geometries)
            for lon, lat in points
        )
        total = len(points)
        fraction = pal / total
        realm_ok = pal >= int(minimum_localities) and fraction >= float(minimum_fraction)
        terrestrial = terrestrial_core(
            class_name=class_name,
            order=order,
            family=family,
            allowed_orders=allowed_orders,
            excluded_coleoptera_families=excluded_coleoptera_families,
        )
        summaries.append(
            SpeciesRealmSummary(
                species=species,
                class_name=class_name,
                order=order,
                family=family,
                total_localities=total,
                palearctic_localities=int(pal),
                palearctic_fraction=float(fraction),
                terrestrial_core=bool(terrestrial),
                realm_eligible=bool(realm_ok),
                biologically_eligible=bool(terrestrial and realm_ok),
            )
        )
    return tuple(summaries)


def apply_gbif_admissibility(
    summaries: Sequence[SpeciesRealmSummary],
    retained_occurrences: Mapping[str, int],
    *,
    minimum_occurrences: int,
) -> tuple[str, ...]:
    eligible = {row.species for row in summaries if row.biologically_eligible}
    missing = sorted(eligible - set(retained_occurrences))
    if missing:
        raise ValueError(f"GBIF admissibility ledger missing biologically eligible species: {missing[:5]}")
    return tuple(
        sorted(
            species
            for species in eligible
            if int(retained_occurrences[species]) >= int(minimum_occurrences)
        )
    )


def deterministic_role_split(
    species: Sequence[str],
    *,
    namespace: str,
) -> tuple[tuple[str, ...], tuple[str, ...]]:
    names = tuple(sorted(map(str, species)))
    if len(set(names)) != len(names):
        raise ValueError("role-split species must be unique")
    keyed = sorted(
        (
            hashlib.sha256(f"{namespace}|{name}".encode("utf-8")).hexdigest(),
            name,
        )
        for name in names
    )
    n_source = len(keyed) // 2
    source = tuple(sorted(name for _, name in keyed[:n_source]))
    target = tuple(sorted(name for _, name in keyed[n_source:]))
    return source, target


def load_geojson(path: Path) -> dict:
    return json.loads(Path(path).read_text(encoding="utf-8"))


__all__ = [
    "SpeciesRealmSummary",
    "apply_gbif_admissibility",
    "deterministic_role_split",
    "geometry_contains_or_boundary",
    "load_geojson",
    "point_in_any_geometry",
    "select_realm_geometries",
    "summarize_species_realm",
    "terrestrial_core",
]
