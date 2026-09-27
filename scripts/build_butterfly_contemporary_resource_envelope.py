#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import json
from collections import defaultdict
from pathlib import Path

from ttf.lepidoptera_host_resource import build_insect_host_footprints


def load_sidecars(interactions: Path, distributions: Path):
    pairs: list[tuple[str, str]] = []
    with interactions.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        required = {"insect_species", "accepted_plant_name_id"}
        if not required <= set(reader.fieldnames or ()):
            raise RuntimeError("insect-host sidecar schema drift")
        for row in reader:
            insect = str(row["insect_species"]).strip()
            host = str(row["accepted_plant_name_id"]).strip()
            if insect and host:
                pairs.append((insect, host))

    units: dict[str, set[str]] = defaultdict(set)
    with distributions.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        required = {"accepted_plant_name_id", "area_code_l3"}
        if not required <= set(reader.fieldnames or ()):
            raise RuntimeError("host-distribution sidecar schema drift")
        for row in reader:
            host = str(row["accepted_plant_name_id"]).strip()
            unit = str(row["area_code_l3"]).strip()
            if host and unit:
                units[host].add(unit)
    return build_insect_host_footprints(pairs, units)


def main() -> int:
    ap = argparse.ArgumentParser()
    panel_group = ap.add_mutually_exclusive_group(required=True)
    panel_group.add_argument("--pilot-json", type=Path)
    panel_group.add_argument("--panel-json", type=Path)
    ap.add_argument("--native-footprints-json", type=Path, required=True)
    ap.add_argument("--insect-host-csv", type=Path, required=True)
    ap.add_argument("--contemporary-distribution-csv", type=Path, required=True)
    ap.add_argument("--output-footprints", type=Path, required=True)
    ap.add_argument("--output-comparison", type=Path, required=True)
    args = ap.parse_args()

    panel_path = args.panel_json if args.panel_json is not None else args.pilot_json
    panel = json.loads(panel_path.read_text(encoding="utf-8"))
    schema = panel.get("schema")
    if schema == "ttf_butterfly_resource_envelope_pilot_v0.1":
        species = tuple(map(str, panel["pilot_species"]))
        if len(species) != 10 or len(set(species)) != 10:
            raise RuntimeError("expected exact ten-species pilot")
    elif schema == "ttf_butterfly_climate_release_independent_panel_v0.1":
        if panel.get("status") != "FROZEN_BEFORE_INDEPENDENT_GBIF_OR_CLIMATE":
            raise RuntimeError("independent panel is not frozen")
        species = tuple(map(str, panel["species"]))
        if len(species) != 32 or len(set(species)) != 32:
            raise RuntimeError("expected exact 32-species independent panel")
    else:
        raise RuntimeError("unexpected butterfly species-panel schema")

    native_payload = json.loads(args.native_footprints_json.read_text(encoding="utf-8"))
    if native_payload.get("schema") != "ttf_butterfly_resource_envelope_s1_footprints_v0.1":
        raise RuntimeError("unexpected native footprint schema")
    native = {
        str(name): frozenset(map(str, units))
        for name, units in native_payload["species"].items()
    }

    contemporary, diagnostics = load_sidecars(
        args.insect_host_csv,
        args.contemporary_distribution_csv,
    )

    rows = []
    for name in species:
        n = native.get(name, frozenset())
        c = contemporary.get(name, frozenset())
        if not n.issubset(c):
            raise RuntimeError(f"native units not subset of contemporary units: {name}")
        rows.append(
            {
                "species": name,
                "native_units": len(n),
                "contemporary_units": len(c),
                "introduced_added_units": len(c - n),
                "contemporary_over_native_ratio": (
                    None if not n else len(c) / len(n)
                ),
                "resolved_host_species": int(
                    diagnostics.get(name, {}).get("resolved_host_species", 0)
                ),
                "hosts_with_contemporary_units": int(
                    diagnostics.get(name, {}).get(
                        "hosts_with_primary_native_units", 0
                    )
                ),
            }
        )

    footprints = {
        "schema": "ttf_butterfly_resource_envelope_contemporary_footprints_v0.1",
        "status": "EXPLORATORY_CONTEMPORARY_HOST_RESOURCE_RECONSTRUCTION",
        "definition": (
            "Union of WCVP-v13 WGSRPD3 units for resolved larval host species "
            "with extinct=0 and location_doubtful=0, retaining both native and "
            "introduced distribution records."
        ),
        "panel_species": list(species),
        "species": {
            name: sorted(contemporary.get(name, frozenset()))
            for name in species
        },
        "claim_boundary": {
            "contemporary_means_wcvp_extant_nondoubtful_not_time_stamped_2010_2026": True,
            "introduced_ranges_included": True,
            "host_database_completeness_not_proven": True,
            "butterfly_occurrences_used": False,
            "genetic_response_used": False,
        },
    }
    comparison = {
        "schema": "ttf_butterfly_resource_envelope_native_contemporary_comparison_v0.1",
        "status": "EXPLORATORY_RESOURCE_ENVELOPE_SENSITIVITY",
        "species": rows,
        "summary": {
            "panel_species": len(species),
            "species_expanded_by_introduced_host_ranges": sum(
                int(row["introduced_added_units"] > 0) for row in rows
            ),
            "total_native_species_units": sum(row["native_units"] for row in rows),
            "total_contemporary_species_units": sum(
                row["contemporary_units"] for row in rows
            ),
            "total_introduced_added_species_units": sum(
                row["introduced_added_units"] for row in rows
            ),
        },
        "interpretation_boundary": footprints["claim_boundary"],
    }

    args.output_footprints.parent.mkdir(parents=True, exist_ok=True)
    args.output_footprints.write_text(
        json.dumps(footprints, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    args.output_comparison.write_text(
        json.dumps(comparison, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    print(json.dumps(comparison, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
