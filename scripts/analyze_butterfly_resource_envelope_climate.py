#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import json
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np
from shapely.geometry import Point, shape
from shapely.strtree import STRtree

from ttf.resource_envelope_climate import (
    climate_mismatch,
    deterministic_interior_points,
    observed_unit_split,
    probability_greater,
)


VARIABLES = ("bio1", "bio7", "bio12", "bio15")


def load_level3(path: Path):
    payload = json.loads(path.read_text(encoding="utf-8"))
    geometries = []
    codes = []
    names = []
    for feature in payload.get("features") or []:
        props = feature.get("properties") or {}
        code = str(props.get("LEVEL3_COD") or "").strip()
        name = str(props.get("LEVEL3_NAM") or "").strip()
        geom = feature.get("geometry")
        if code and geom:
            geometries.append(shape(geom))
            codes.append(code)
            names.append(name)
    if not geometries or len(codes) != len(set(codes)):
        raise RuntimeError("WGSRPD3 geometry/code drift")
    return geometries, tuple(codes), tuple(names)


def point_mapper(geometries, codes):
    tree = STRtree(geometries)
    index_by_id = {id(geom): i for i, geom in enumerate(geometries)}

    def map_point(lon: float, lat: float):
        point = Point(float(lon), float(lat))
        matches = []
        for candidate in tree.query(point):
            if hasattr(candidate, "geom_type"):
                idx = index_by_id[id(candidate)]
            else:
                idx = int(candidate)
            if geometries[idx].covers(point):
                matches.append(codes[idx])
        if not matches:
            return None
        return sorted(set(matches))[0]

    return map_point


def load_unit_table(path: Path):
    rows = []
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        required = {
            "species",
            "wgsrpd3_code",
            "host_available",
            "butterfly_observed",
            "butterfly_record_count",
        }
        if not required <= set(reader.fieldnames or ()):
            raise RuntimeError("WGSRPD3 unit table schema drift")
        for row in reader:
            rows.append(
                {
                    "species": str(row["species"]).strip(),
                    "unit": str(row["wgsrpd3_code"]).strip(),
                    "host_available": int(row["host_available"]),
                    "butterfly_observed": int(row["butterfly_observed"]),
                    "butterfly_record_count": int(row["butterfly_record_count"]),
                }
            )
    return rows


def sample_matrix(datasets, coords):
    columns = []
    valid = np.ones(len(coords), dtype=bool)
    for dataset in datasets:
        values = np.asarray(
            [float(v[0]) for v in dataset.sample(coords)],
            dtype=float,
        )
        if dataset.nodata is not None:
            valid &= ~np.isclose(values, float(dataset.nodata), rtol=0.0, atol=0.0)
        valid &= np.isfinite(values)
        columns.append(values)
    return np.column_stack(columns), valid


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--occurrences-csv", type=Path, required=True)
    ap.add_argument("--unit-table", type=Path, required=True)
    ap.add_argument("--contemporary-footprints-json", type=Path, required=True)
    ap.add_argument("--gate-rule-json", type=Path, required=True)
    ap.add_argument("--quality-gate-json", type=Path, required=True)
    ap.add_argument("--level3-geojson", type=Path, required=True)
    ap.add_argument("--raster", type=Path, action="append", required=True)
    ap.add_argument("--unit-sample-points", type=int, default=16)
    ap.add_argument("--output-unit-climate", type=Path, required=True)
    ap.add_argument("--output-summary", type=Path, required=True)
    args = ap.parse_args()

    if len(args.raster) != 4:
        raise RuntimeError("exactly four CHELSA rasters are required in bio1,bio7,bio12,bio15 order")

    gate_rule = json.loads(args.gate_rule_json.read_text(encoding="utf-8"))
    quality = json.loads(args.quality_gate_json.read_text(encoding="utf-8"))
    gate_schema = str(gate_rule.get("schema") or "")
    quality_schema = str(quality.get("schema") or "")

    if gate_schema == "ttf_butterfly_resource_envelope_climate_pilot_gate_v0.1":
        if gate_rule.get("status") != "FROZEN_EXPLORATORY_QUALITY_GATE_BEFORE_ANY_CLIMATE_RESULT":
            raise RuntimeError("climate pilot gate rule is not frozen")
        if quality_schema != "ttf_butterfly_resource_envelope_preclimate_quality_gate_result_v0.1":
            raise RuntimeError("unexpected pilot preclimate quality-gate schema")
        if quality.get("status") != "PASS_TO_EXPLORATORY_CLIMATE_PILOT":
            raise RuntimeError("preclimate quality gate did not authorize climate pilot")
        minimum_species = int(
            gate_rule["pilot_level_gate"]["minimum_species_passing_quality_gate"]
        )
        q = gate_rule["species_quality_gate"]
        effort_threshold = int(
            q["sampling_effort_identifiability"]["other_pilot_record_threshold"]
        )
        minimum_training_records = int(
            q["climate_training_floor"]["minimum_valid_training_occurrence_records"]
        )
        minimum_eval_units = 2
        minimum_never_units = 2
        analysis_scope = "response_blind_pilot"
    elif gate_schema in {
        "ttf_butterfly_climate_release_independent_test_v0.1",
        "ttf_butterfly_climate_release_independent_test_v0.2",
    }:
        expected_status = {
            "ttf_butterfly_climate_release_independent_test_v0.1": (
                "FROZEN_PILOT_DERIVED_HYPOTHESIS_BEFORE_INDEPENDENT_GBIF_OR_CLIMATE"
            ),
            "ttf_butterfly_climate_release_independent_test_v0.2": (
                "FROZEN_RESPONSE_BLIND_PRIMARY_INFERENCE_CORRECTION_BEFORE_INDEPENDENT_PRECLIMATE_OR_CLIMATE_RESULT"
            ),
        }[gate_schema]
        if gate_rule.get("status") != expected_status:
            raise RuntimeError("independent climate-release protocol is not frozen")
        if quality_schema != "ttf_butterfly_climate_release_preclimate_gate_v0.1":
            raise RuntimeError("unexpected independent preclimate quality-gate schema")
        if quality.get("status") != "PASS_TO_INDEPENDENT_CLIMATE_CROSSFIT":
            raise RuntimeError("independent preclimate gate did not authorize cross-fit")
        minimum_species = int(
            gate_rule["preclimate_quality_gate"][
                "minimum_species_passing_preclimate_gate"
            ]
        )
        q = gate_rule["preclimate_quality_gate"]["thresholds"]
        effort_threshold = int(q["other_panel_record_effort_threshold"])
        climate_rule = gate_rule["climate_crossfit"]
        minimum_training_records = int(
            climate_rule["minimum_training_occurrence_records"]
        )
        minimum_eval_units = int(
            climate_rule["minimum_effort_supported_eval_observed_units"]
        )
        minimum_never_units = int(
            climate_rule["minimum_effort_supported_never_observed_units"]
        )
        analysis_scope = "independent_pilot_derived_test"
    else:
        raise RuntimeError("unexpected climate analysis protocol schema")

    species = sorted(map(str, quality.get("qualified_species", [])))
    if len(species) < minimum_species:
        raise RuntimeError("qualified species fell below frozen minimum")
    if (
        effort_threshold < 0
        or minimum_training_records < 2
        or minimum_eval_units < 1
        or minimum_never_units < 1
    ):
        raise RuntimeError("invalid frozen climate quality threshold")

    contemporary_payload = json.loads(
        args.contemporary_footprints_json.read_text(encoding="utf-8")
    )
    if contemporary_payload.get("schema") != (
        "ttf_butterfly_resource_envelope_contemporary_footprints_v0.1"
    ):
        raise RuntimeError("unexpected contemporary host-footprint schema")
    contemporary = {
        str(name): frozenset(map(str, units))
        for name, units in contemporary_payload.get("species", {}).items()
    }
    missing_contemporary = [name for name in species if name not in contemporary]
    if missing_contemporary:
        raise RuntimeError(
            "qualified species missing contemporary host footprint: "
            + ", ".join(missing_contemporary)
        )

    unit_rows = load_unit_table(args.unit_table)
    geometries, codes, code_names = load_level3(args.level3_geojson)
    code_to_geom = dict(zip(codes, geometries))
    code_to_name = dict(zip(codes, code_names))
    map_point = point_mapper(geometries, codes)

    observed_by_species = defaultdict(set)
    records_by_unit_species = defaultdict(Counter)
    for row in unit_rows:
        name = str(row["species"])
        code = str(row["unit"])
        if int(row["butterfly_observed"]):
            observed_by_species[name].add(code)
        count = int(row["butterfly_record_count"])
        if count:
            records_by_unit_species[code][name] += count

    host_by_species = {
        name: set(contemporary[name])
        for name in species
    }
    required_units = sorted(set().union(*(host_by_species[name] for name in species)))
    missing_geometry = [code for code in required_units if code not in code_to_geom]
    if missing_geometry:
        raise RuntimeError(
            "host units missing WGSRPD3 geometry: " + ", ".join(missing_geometry[:10])
        )

    import rasterio

    datasets = [rasterio.open(path) for path in args.raster]
    try:
        for dataset in datasets:
            if dataset.crs is None or not dataset.crs.is_geographic:
                raise RuntimeError(f"CHELSA raster is not geographic: {dataset.name}")

        unit_climate = {}
        unit_point_n = {}
        for code in required_units:
            coords = deterministic_interior_points(
                code_to_geom[code],
                maximum_points=int(args.unit_sample_points),
            )
            if not coords:
                continue
            matrix, valid = sample_matrix(datasets, coords)
            matrix = matrix[valid]
            if len(matrix) == 0:
                continue
            unit_climate[code] = np.median(matrix, axis=0)
            unit_point_n[code] = int(len(matrix))

        occurrence_rows = []
        with args.occurrences_csv.open(newline="", encoding="utf-8") as handle:
            reader = csv.DictReader(handle)
            required = {"species", "gbif_key", "latitude", "longitude"}
            if not required <= set(reader.fieldnames or ()):
                raise RuntimeError("occurrence CSV schema drift")
            for row in reader:
                name = str(row["species"]).strip()
                if name not in species:
                    continue
                lon = float(row["longitude"])
                lat = float(row["latitude"])
                code = map_point(lon, lat)
                if code is None:
                    continue
                occurrence_rows.append(
                    {
                        "species": name,
                        "gbif_key": int(row["gbif_key"]),
                        "unit": code,
                        "lon": lon,
                        "lat": lat,
                    }
                )

        occurrence_coords = [(r["lon"], r["lat"]) for r in occurrence_rows]
        occurrence_climate, occurrence_valid = sample_matrix(datasets, occurrence_coords)
        valid_occurrences = []
        for row, vector, valid in zip(
            occurrence_rows,
            occurrence_climate,
            occurrence_valid,
        ):
            if valid:
                valid_occurrences.append((row, np.asarray(vector, dtype=float)))

        occ_by_species_unit = defaultdict(list)
        for row, vector in valid_occurrences:
            occ_by_species_unit[(str(row["species"]), str(row["unit"]))].append(vector)

        result_rows = []
        summaries = []
        for name in species:
            host_units = host_by_species[name]
            observed = observed_by_species[name] & host_units
            train_units, eval_units = observed_unit_split(name, observed)
            train_vectors = []
            for code in train_units:
                train_vectors.extend(occ_by_species_unit.get((name, code), []))
            train = np.asarray(train_vectors, dtype=float)

            never_observed = host_units - observed
            eval_host_observed = set(eval_units)

            eligible_codes = sorted(
                code
                for code in (never_observed | eval_host_observed)
                if code in unit_climate
                and sum(
                    count
                    for other, count in records_by_unit_species.get(code, {}).items()
                    if other != name
                ) >= int(effort_threshold)
            )

            mismatch_by_code = {}
            if len(train) >= 2 and eligible_codes:
                matrix = np.vstack([unit_climate[code] for code in eligible_codes])
                try:
                    mismatch = climate_mismatch(matrix, train)
                except ValueError:
                    mismatch = np.full(len(eligible_codes), np.nan)
                mismatch_by_code = dict(zip(eligible_codes, map(float, mismatch)))

            occupied_mismatch = [
                mismatch_by_code[code]
                for code in sorted(eval_host_observed)
                if code in mismatch_by_code and np.isfinite(mismatch_by_code[code])
            ]
            unoccupied_mismatch = [
                mismatch_by_code[code]
                for code in sorted(never_observed)
                if code in mismatch_by_code and np.isfinite(mismatch_by_code[code])
            ]

            for code in sorted(host_units):
                other_effort = sum(
                    count
                    for other, count in records_by_unit_species.get(code, {}).items()
                    if other != name
                )
                result_rows.append(
                    {
                        "species": name,
                        "wgsrpd3_code": code,
                        "wgsrpd3_name": code_to_name.get(code, ""),
                        "observed_unit_split": (
                            "train_observed"
                            if code in train_units
                            else "eval_observed"
                            if code in eval_host_observed
                            else "never_observed"
                        ),
                        "other_pilot_record_effort": other_effort,
                        "effort_supported": int(
                            other_effort >= int(effort_threshold)
                        ),
                        "climate_sample_points": unit_point_n.get(code, 0),
                        "bio1": (
                            None if code not in unit_climate else float(unit_climate[code][0])
                        ),
                        "bio7": (
                            None if code not in unit_climate else float(unit_climate[code][1])
                        ),
                        "bio12": (
                            None if code not in unit_climate else float(unit_climate[code][2])
                        ),
                        "bio15": (
                            None if code not in unit_climate else float(unit_climate[code][3])
                        ),
                        "climate_mismatch_to_train_niche": mismatch_by_code.get(code),
                    }
                )

            delta = None
            if occupied_mismatch and unoccupied_mismatch:
                delta = float(np.median(unoccupied_mismatch) - np.median(occupied_mismatch))
            summaries.append(
                {
                    "species": name,
                    "host_units": len(host_units),
                    "observed_units": len(observed),
                    "train_observed_units": len(train_units),
                    "eval_observed_host_units": len(eval_host_observed),
                    "never_observed_host_units": len(never_observed),
                    "training_occurrence_records_with_climate": int(len(train)),
                    "effort_threshold_other_pilot_records": int(effort_threshold),
                    "effort_supported_eval_observed_units": len(occupied_mismatch),
                    "effort_supported_never_observed_units": len(unoccupied_mismatch),
                    "median_mismatch_eval_observed": (
                        None
                        if not occupied_mismatch
                        else float(np.median(occupied_mismatch))
                    ),
                    "median_mismatch_never_observed": (
                        None
                        if not unoccupied_mismatch
                        else float(np.median(unoccupied_mismatch))
                    ),
                    "mismatch_delta_never_minus_eval_observed": delta,
                    "probability_never_observed_mismatch_greater": probability_greater(
                        unoccupied_mismatch,
                        occupied_mismatch,
                    ),
                    "climate_crossfit_informative": (
                        len(train) >= minimum_training_records
                        and len(occupied_mismatch) >= minimum_eval_units
                        and len(unoccupied_mismatch) >= minimum_never_units
                    ),
                }
            )
    finally:
        for dataset in datasets:
            dataset.close()

    args.output_unit_climate.parent.mkdir(parents=True, exist_ok=True)
    fields = [
        "species",
        "wgsrpd3_code",
        "wgsrpd3_name",
        "observed_unit_split",
        "other_pilot_record_effort",
        "effort_supported",
        "climate_sample_points",
        "bio1",
        "bio7",
        "bio12",
        "bio15",
        "climate_mismatch_to_train_niche",
    ]
    with args.output_unit_climate.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(result_rows)

    payload = {
        "schema": "ttf_butterfly_resource_envelope_climate_crossfit_v0.1",
        "status": "EXPLORATORY_CLIMATE_CROSSFIT_DIAGNOSTIC",
        "variables": list(VARIABLES),
        "analysis_program_schema": gate_schema,
        "analysis_scope": analysis_scope,
        "quality_qualified_species": list(species),
        "quality_gate_minimum_species": minimum_species,
        "frozen_other_pilot_effort_threshold": effort_threshold,
        "frozen_minimum_training_occurrence_records": minimum_training_records,
        "frozen_minimum_effort_supported_eval_observed_units": minimum_eval_units,
        "frozen_minimum_effort_supported_never_observed_units": minimum_never_units,
        "host_envelope": "WCVP-v13 extant non-doubtful contemporary host-resource WGSRPD3 units including introduced ranges",
        "unit_climate": (
            "Median CHELSA V2.1 value over up to a fixed number of deterministic "
            "interior points per WGSRPD3 polygon/component."
        ),
        "crossfit": (
            "Observed WGSRPD3 units inside the contemporary host-resource envelope "
            "are deterministically hash-split by species. Climate niche center/scale "
            "is estimated only from GBIF occurrences falling in train-observed "
            "within-envelope units. Comparison uses held-out eval-observed host units "
            "versus contemporary host units never observed for that butterfly."
        ),
        "species": summaries,
        "informative_species": sum(
            bool(row["climate_crossfit_informative"]) for row in summaries
        ),
        "claim_boundary": {
            "exploratory_only": True,
            "never_observed_is_not_true_absence": True,
            "other_pilot_effort_is_only_a_sampling_proxy": True,
            "host_database_completeness_not_proven": True,
            "only_preclimate_quality_qualified_species_analyzed": True,
            "contemporary_host_envelope_used": True,
            "crossfit_is_species_level_effect_estimation_not_the_primary_between_species_test": True,
        },
    }
    args.output_summary.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    print(json.dumps(payload, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
