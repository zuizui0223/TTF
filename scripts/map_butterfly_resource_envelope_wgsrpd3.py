#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import json
from collections import Counter, defaultdict
from pathlib import Path

from shapely.geometry import Point, shape
from shapely.strtree import STRtree


def load_species_panel(path: Path) -> tuple[str, ...]:
    payload = json.loads(path.read_text(encoding="utf-8"))
    schema = payload.get("schema")
    if schema == "ttf_butterfly_resource_envelope_pilot_v0.1":
        names = tuple(map(str, payload.get("pilot_species", [])))
        if len(names) != 10 or len(set(names)) != 10:
            raise RuntimeError("expected exact ten-species pilot")
        return names
    if schema == "ttf_butterfly_climate_release_independent_panel_v0.1":
        if payload.get("status") != "FROZEN_BEFORE_INDEPENDENT_GBIF_OR_CLIMATE":
            raise RuntimeError("independent panel is not frozen")
        names = tuple(map(str, payload.get("species", [])))
        if len(names) != 32 or len(set(names)) != 32:
            raise RuntimeError("expected exact 32-species independent panel")
        return names
    raise RuntimeError("unexpected butterfly species-panel schema")


def load_pilot(path: Path) -> tuple[str, ...]:
    return load_species_panel(path)


def load_host_footprints(path: Path) -> dict[str, frozenset[str]]:
    payload = json.loads(path.read_text(encoding="utf-8"))
    if payload.get("schema") != "ttf_butterfly_resource_envelope_s1_footprints_v0.1":
        raise RuntimeError("unexpected host-footprint schema")
    return {
        str(name): frozenset(map(str, units))
        for name, units in payload.get("species", {}).items()
    }


def load_level3(path: Path):
    payload = json.loads(path.read_text(encoding="utf-8"))
    features = payload.get("features") or []
    geometries = []
    codes = []
    names = []
    for feature in features:
        props = feature.get("properties") or {}
        code = str(props.get("LEVEL3_COD") or "").strip()
        name = str(props.get("LEVEL3_NAM") or "").strip()
        geom_payload = feature.get("geometry")
        if not code or not geom_payload:
            continue
        geometries.append(shape(geom_payload))
        codes.append(code)
        names.append(name)
    if not geometries or len(set(codes)) != len(codes):
        raise RuntimeError("WGSRPD level-3 geometry/code drift")
    return geometries, tuple(codes), tuple(names)


def tree_candidate_indices(tree: STRtree, geometries, point: Point):
    candidates = tree.query(point)
    index_by_id = {id(geom): i for i, geom in enumerate(geometries)}
    out = []
    for candidate in candidates:
        if hasattr(candidate, "geom_type"):
            out.append(index_by_id[id(candidate)])
        else:
            out.append(int(candidate))
    return out


def map_point(tree, geometries, codes, lon: float, lat: float):
    point = Point(float(lon), float(lat))
    matches = []
    for idx in tree_candidate_indices(tree, geometries, point):
        geom = geometries[idx]
        if geom.covers(point):
            matches.append(codes[idx])
    if not matches:
        return None, 0
    matches = sorted(set(matches))
    return matches[0], len(matches)


def main() -> int:
    ap = argparse.ArgumentParser()
    panel_group = ap.add_mutually_exclusive_group(required=True)
    panel_group.add_argument("--pilot-json", type=Path)
    panel_group.add_argument("--panel-json", type=Path)
    ap.add_argument("--host-footprints-json", type=Path, required=True)
    ap.add_argument("--occurrences-csv", type=Path, required=True)
    ap.add_argument("--level3-geojson", type=Path, required=True)
    ap.add_argument("--output-unit-table", type=Path, required=True)
    ap.add_argument("--output-summary", type=Path, required=True)
    args = ap.parse_args()

    panel_path = args.panel_json if args.panel_json is not None else args.pilot_json
    pilot = load_species_panel(panel_path)
    footprints = load_host_footprints(args.host_footprints_json)
    missing_footprints = [name for name in pilot if name not in footprints]
    if missing_footprints:
        raise RuntimeError(
            "panel species missing host footprints: " + ", ".join(missing_footprints)
        )

    geometries, codes, code_names = load_level3(args.level3_geojson)
    tree = STRtree(geometries)
    code_to_name = dict(zip(codes, code_names))

    record_counts: dict[str, Counter[str]] = defaultdict(Counter)
    raw_records = Counter()
    unmapped_records = Counter()
    multiply_mapped_records = Counter()
    with args.occurrences_csv.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        required = {"species", "gbif_key", "latitude", "longitude"}
        if not required <= set(reader.fieldnames or ()):
            raise RuntimeError("occurrence CSV schema drift")
        for row in reader:
            species = str(row["species"]).strip()
            if species not in pilot:
                continue
            raw_records[species] += 1
            code, n_matches = map_point(
                tree,
                geometries,
                codes,
                float(row["longitude"]),
                float(row["latitude"]),
            )
            if code is None:
                unmapped_records[species] += 1
                continue
            if n_matches > 1:
                multiply_mapped_records[species] += 1
            record_counts[species][code] += 1

    unit_rows = []
    species_summaries = []
    informative = 0
    for species in pilot:
        host_units = set(footprints[species])
        observed_units = set(record_counts[species])
        within = observed_units & host_units
        outside = observed_units - host_units
        unobserved_host = host_units - observed_units
        resolution_informative = (
            len(host_units) >= 5
            and len(within) >= 2
            and len(unobserved_host) >= 2
        )
        informative += int(resolution_informative)

        for code in sorted(host_units | observed_units):
            unit_rows.append(
                {
                    "species": species,
                    "wgsrpd3_code": code,
                    "wgsrpd3_name": code_to_name.get(code, ""),
                    "host_available": int(code in host_units),
                    "butterfly_observed": int(code in observed_units),
                    "butterfly_record_count": int(record_counts[species].get(code, 0)),
                }
            )

        species_summaries.append(
            {
                "species": species,
                "host_units": len(host_units),
                "butterfly_observed_units": len(observed_units),
                "observed_within_host_units": len(within),
                "observed_outside_host_units": len(outside),
                "host_units_without_sampled_butterfly": len(unobserved_host),
                "raw_resource_fill_fraction": (
                    None if not host_units else len(within) / len(host_units)
                ),
                "raw_occurrence_records": int(raw_records[species]),
                "records_unmapped_to_wgsrpd3": int(unmapped_records[species]),
                "records_on_multiple_polygon_boundaries": int(
                    multiply_mapped_records[species]
                ),
                "resolution_informative": resolution_informative,
            }
        )

    args.output_unit_table.parent.mkdir(parents=True, exist_ok=True)
    fields = [
        "species",
        "wgsrpd3_code",
        "wgsrpd3_name",
        "host_available",
        "butterfly_observed",
        "butterfly_record_count",
    ]
    with args.output_unit_table.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(unit_rows)

    summary = {
        "schema": "ttf_butterfly_resource_envelope_wgsrpd3_pilot_v0.1",
        "status": "EXPLORATORY_RESOLUTION_DIAGNOSTIC",
        "pilot_species": list(pilot),
        "species": species_summaries,
        "resolution_summary": {
            "informative_species": informative,
            "pilot_species": len(pilot),
            "informative_definition": (
                "host_units>=5 AND observed_within_host_units>=2 AND "
                "host_units_without_sampled_butterfly>=2"
            ),
        },
        "interpretation_boundary": {
            "raw_resource_fill_fraction_is_descriptive_only": True,
            "unobserved_unit_is_not_true_absence": True,
            "no_sampling_effort_correction_yet": True,
            "no_climate_model_yet": True,
            "no_confirmatory_test": True,
        },
    }
    args.output_summary.parent.mkdir(parents=True, exist_ok=True)
    args.output_summary.write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    print(json.dumps(summary["resolution_summary"], indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
