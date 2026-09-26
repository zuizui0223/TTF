#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import json
from collections import defaultdict
from pathlib import Path


THRESHOLDS = (1, 5, 10, 20)


def load_rows(path: Path) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
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
                    "species": str(row["species"]),
                    "unit": str(row["wgsrpd3_code"]),
                    "host_available": int(row["host_available"]),
                    "butterfly_observed": int(row["butterfly_observed"]),
                    "butterfly_record_count": int(row["butterfly_record_count"]),
                }
            )
    return rows


def diagnose(
    rows: list[dict[str, object]],
    *,
    host_footprints: dict[str, frozenset[str]] | None = None,
) -> tuple[list[dict[str, object]], dict]:
    row_species = sorted({str(row["species"]) for row in rows})
    species = (
        sorted(host_footprints)
        if host_footprints is not None
        else row_species
    )
    records_by_unit_species: dict[str, dict[str, int]] = defaultdict(dict)
    observed_by_species: dict[str, set[str]] = defaultdict(set)
    native_host_by_species: dict[str, set[str]] = defaultdict(set)
    for row in rows:
        target = str(row["species"])
        unit = str(row["unit"])
        records_by_unit_species[unit][target] = int(row["butterfly_record_count"])
        if int(row["butterfly_observed"]):
            observed_by_species[target].add(unit)
        if int(row["host_available"]):
            native_host_by_species[target].add(unit)

    out_rows: list[dict[str, object]] = []
    for target in species:
        host_units = (
            set(host_footprints.get(target, frozenset()))
            if host_footprints is not None
            else set(native_host_by_species.get(target, set()))
        )
        for threshold in THRESHOLDS:
            supported = []
            for unit in sorted(host_units):
                other_records = sum(
                    count
                    for sp, count in records_by_unit_species.get(unit, {}).items()
                    if sp != target
                )
                if other_records >= threshold:
                    supported.append((unit, other_records))

            occupied = sum(
                int(unit in observed_by_species.get(target, set()))
                for unit, _ in supported
            )
            unoccupied = len(supported) - occupied
            fill = None if not supported else occupied / len(supported)
            informative = (
                len(supported) >= 5
                and occupied >= 2
                and unoccupied >= 2
            )
            out_rows.append(
                {
                    "species": target,
                    "other_pilot_record_threshold": threshold,
                    "host_units_total": len(host_units),
                    "effort_supported_host_units": len(supported),
                    "occupied_effort_supported_host_units": occupied,
                    "unoccupied_effort_supported_host_units": unoccupied,
                    "descriptive_fill_fraction": fill,
                    "effort_sensitivity_informative": informative,
                }
            )

    by_threshold = {}
    for threshold in THRESHOLDS:
        subset = [
            row
            for row in out_rows
            if int(row["other_pilot_record_threshold"]) == threshold
        ]
        by_threshold[str(threshold)] = {
            "species_with_effort_supported_host_units": sum(
                int(row["effort_supported_host_units"]) > 0 for row in subset
            ),
            "informative_species": sum(
                bool(row["effort_sensitivity_informative"]) for row in subset
            ),
        }

    summary = {
        "schema": "ttf_butterfly_resource_envelope_sampling_effort_sensitivity_v0.1",
        "status": "EXPLORATORY_SAMPLING_EFFORT_DIAGNOSTIC",
        "host_envelope": (
            "contemporary_extant_nondoubtful_including_introduced"
            if host_footprints is not None
            else "native_only_from_unit_table"
        ),
        "effort_proxy": (
            "For each target species and WGSRPD3 unit, the total GBIF record count "
            "of the other nine frozen pilot species in that unit."
        ),
        "thresholds": list(THRESHOLDS),
        "informative_definition": (
            "Among host-available units meeting the other-pilot-record threshold: "
            "at least 5 supported units, at least 2 target-observed units, and at "
            "least 2 target-unobserved units."
        ),
        "by_threshold": by_threshold,
        "interpretation_boundary": {
            "proxy_is_not_global_butterfly_sampling_effort": True,
            "unobserved_unit_is_not_true_absence": True,
            "descriptive_fill_fraction_is_not_occupancy_probability": True,
            "no_climate_model_yet": True,
            "no_confirmatory_test": True,
        },
    }
    return out_rows, summary


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--unit-table", type=Path, required=True)
    ap.add_argument("--host-footprints-json", type=Path)
    ap.add_argument("--output-csv", type=Path, required=True)
    ap.add_argument("--output-json", type=Path, required=True)
    args = ap.parse_args()

    host_footprints = None
    if args.host_footprints_json is not None:
        payload = json.loads(args.host_footprints_json.read_text(encoding="utf-8"))
        if payload.get("schema") != (
            "ttf_butterfly_resource_envelope_contemporary_footprints_v0.1"
        ):
            raise RuntimeError("unexpected contemporary host-footprint schema")
        host_footprints = {
            str(name): frozenset(map(str, units))
            for name, units in payload.get("species", {}).items()
        }
    rows, summary = diagnose(
        load_rows(args.unit_table),
        host_footprints=host_footprints,
    )
    args.output_csv.parent.mkdir(parents=True, exist_ok=True)
    fields = [
        "species",
        "other_pilot_record_threshold",
        "host_units_total",
        "effort_supported_host_units",
        "occupied_effort_supported_host_units",
        "unoccupied_effort_supported_host_units",
        "descriptive_fill_fraction",
        "effort_sensitivity_informative",
    ]
    with args.output_csv.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)
    args.output_json.write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    print(json.dumps(summary, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
