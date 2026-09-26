#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--gate-json", type=Path, required=True)
    ap.add_argument("--descriptors-csv", type=Path, required=True)
    ap.add_argument("--overlap-json", type=Path, required=True)
    ap.add_argument("--effort-csv", type=Path, required=True)
    ap.add_argument("--output-json", type=Path, required=True)
    args = ap.parse_args()

    gate = json.loads(args.gate_json.read_text(encoding="utf-8"))
    if gate.get("schema") != "ttf_butterfly_resource_envelope_climate_pilot_gate_v0.1":
        raise RuntimeError("unexpected climate pilot gate schema")
    q = gate["species_quality_gate"]
    effort_threshold = int(
        q["sampling_effort_identifiability"]["other_pilot_record_threshold"]
    )
    minimum_units = int(
        q["sampling_effort_identifiability"]["minimum_effort_supported_host_units"]
    )
    minimum_observed = int(
        q["sampling_effort_identifiability"][
            "minimum_butterfly_observed_effort_supported_host_units"
        ]
    )
    minimum_unobserved = int(
        q["sampling_effort_identifiability"][
            "minimum_butterfly_unobserved_effort_supported_host_units"
        ]
    )

    descriptors = {}
    with args.descriptors_csv.open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            descriptors[str(row["species"])] = row

    overlap_payload = json.loads(args.overlap_json.read_text(encoding="utf-8"))
    if overlap_payload.get("schema") != (
        "ttf_butterfly_resource_envelope_native_contemporary_occurrence_overlap_v0.1"
    ):
        raise RuntimeError("unexpected overlap diagnostic schema")
    overlap = {str(row["species"]): row for row in overlap_payload["species"]}

    effort = {}
    with args.effort_csv.open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            if int(row["other_pilot_record_threshold"]) == effort_threshold:
                effort[str(row["species"])] = row

    species = sorted(set(overlap) | set(effort))
    rows = []
    qualified = []
    for name in species:
        descriptor = descriptors.get(name)
        overlap_row = overlap.get(name)
        effort_row = effort.get(name)
        reasons = []
        if descriptor is None:
            reasons.append("MISSING_RESOURCE_DESCRIPTOR")
        if overlap_row is None:
            reasons.append("MISSING_CONTEMPORARY_OVERLAP")
        if effort_row is None:
            reasons.append("MISSING_EFFORT_DIAGNOSTIC")

        if descriptor is not None:
            host_families = float(descriptor["host_family_count"])
            resolved_hosts = int(descriptor["resolved_host_species"])
            host_ok = resolved_hosts >= host_families
            if not host_ok:
                reasons.append("HOST_TAXONOMY_LOWER_BOUND_FAIL")
        else:
            host_families = None
            resolved_hosts = None
            host_ok = False

        if overlap_row is not None:
            observed_units = int(overlap_row["butterfly_observed_units"])
            outside = int(overlap_row["observed_outside_contemporary"])
            outside_fraction = None if observed_units == 0 else outside / observed_units
            envelope_ok = (
                outside_fraction is not None
                and outside_fraction
                <= float(
                    q["contemporary_envelope_overlap"][
                        "maximum_outside_contemporary_fraction"
                    ]
                )
            )
            if not envelope_ok:
                reasons.append("CONTEMPORARY_ENVELOPE_OVERLAP_FAIL")
        else:
            observed_units = 0
            outside = 0
            outside_fraction = None
            envelope_ok = False

        if effort_row is not None:
            supported = int(effort_row["effort_supported_host_units"])
            observed_supported = int(
                effort_row["occupied_effort_supported_host_units"]
            )
            unobserved_supported = int(
                effort_row["unoccupied_effort_supported_host_units"]
            )
            effort_ok = (
                supported >= minimum_units
                and observed_supported >= minimum_observed
                and unobserved_supported >= minimum_unobserved
            )
            if not effort_ok:
                reasons.append("SAMPLING_EFFORT_IDENTIFIABILITY_FAIL")
        else:
            supported = observed_supported = unobserved_supported = 0
            effort_ok = False

        passed = not reasons
        if passed:
            qualified.append(name)
        rows.append(
            {
                "species": name,
                "passed_preclimate_quality_gate": passed,
                "host_family_count": host_families,
                "resolved_host_species": resolved_hosts,
                "host_taxonomy_lower_bound_pass": host_ok,
                "butterfly_observed_units": observed_units,
                "outside_contemporary_units": outside,
                "outside_contemporary_fraction": outside_fraction,
                "contemporary_envelope_overlap_pass": envelope_ok,
                "effort_threshold": effort_threshold,
                "effort_supported_host_units": supported,
                "effort_supported_observed_units": observed_supported,
                "effort_supported_unobserved_units": unobserved_supported,
                "sampling_effort_identifiability_pass": effort_ok,
                "failure_reasons": reasons,
            }
        )

    minimum_species = int(gate["pilot_level_gate"]["minimum_species_passing_quality_gate"])
    payload = {
        "schema": "ttf_butterfly_resource_envelope_preclimate_quality_gate_result_v0.1",
        "status": (
            "PASS_TO_EXPLORATORY_CLIMATE_PILOT"
            if len(qualified) >= minimum_species
            else str(gate["pilot_level_gate"]["if_below"])
        ),
        "qualified_species": qualified,
        "qualified_species_count": len(qualified),
        "minimum_species_required": minimum_species,
        "species": rows,
        "climate_result_used": False,
        "genetic_response_used": False,
    }
    args.output_json.parent.mkdir(parents=True, exist_ok=True)
    args.output_json.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    print(json.dumps(payload, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
