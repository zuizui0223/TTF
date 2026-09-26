#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path


PROTOCOL_SCHEMAS = {
    "ttf_butterfly_climate_release_independent_test_v0.1": (
        "FROZEN_PILOT_DERIVED_HYPOTHESIS_BEFORE_INDEPENDENT_GBIF_OR_CLIMATE"
    ),
    "ttf_butterfly_climate_release_independent_test_v0.2": (
        "FROZEN_RESPONSE_BLIND_PRIMARY_INFERENCE_CORRECTION_BEFORE_INDEPENDENT_PRECLIMATE_OR_CLIMATE_RESULT"
    ),
}
PANEL_SCHEMA = "ttf_butterfly_climate_release_independent_panel_v0.1"
OVERLAP_SCHEMA = (
    "ttf_butterfly_resource_envelope_native_contemporary_occurrence_overlap_v0.1"
)
LEDGER_MATRIX_SCHEMA = "ttf_butterfly_resource_envelope_occurrence_matrix_v0.1"


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--protocol-json", type=Path, required=True)
    ap.add_argument("--panel-json", type=Path, required=True)
    ap.add_argument("--descriptors-csv", type=Path, required=True)
    ap.add_argument("--overlap-json", type=Path, required=True)
    ap.add_argument("--effort-csv", type=Path, required=True)
    ap.add_argument("--occurrence-ledgers-json", type=Path, required=True)
    ap.add_argument("--output-json", type=Path, required=True)
    args = ap.parse_args()

    protocol = json.loads(args.protocol_json.read_text(encoding="utf-8"))
    protocol_schema = str(protocol.get("schema") or "")
    expected_status = PROTOCOL_SCHEMAS.get(protocol_schema)
    if expected_status is None:
        raise RuntimeError("unexpected independent-test protocol schema")
    if protocol.get("status") != expected_status:
        raise RuntimeError("independent-test protocol is not frozen")

    panel = json.loads(args.panel_json.read_text(encoding="utf-8"))
    if panel.get("schema") != PANEL_SCHEMA:
        raise RuntimeError("unexpected independent-panel schema")
    if panel.get("status") != "FROZEN_BEFORE_INDEPENDENT_GBIF_OR_CLIMATE":
        raise RuntimeError("independent panel is not frozen")
    species = tuple(map(str, panel.get("species", [])))
    if len(species) != 32 or len(set(species)) != 32:
        raise RuntimeError("expected exact 32-species independent panel")

    q = protocol["preclimate_quality_gate"]
    thresholds = q["thresholds"]
    maximum_outside = float(
        thresholds["maximum_outside_contemporary_fraction"]
    )
    effort_threshold = int(
        thresholds["other_panel_record_effort_threshold"]
    )
    minimum_units = int(
        thresholds["minimum_effort_supported_host_units"]
    )
    minimum_observed = int(
        thresholds["minimum_effort_supported_observed_units"]
    )
    minimum_unobserved = int(
        thresholds["minimum_effort_supported_unobserved_units"]
    )
    minimum_species = int(q["minimum_species_passing_preclimate_gate"])

    descriptors = {}
    with args.descriptors_csv.open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            descriptors[str(row["species"])] = row

    overlap_payload = json.loads(args.overlap_json.read_text(encoding="utf-8"))
    if overlap_payload.get("schema") != OVERLAP_SCHEMA:
        raise RuntimeError("unexpected overlap diagnostic schema")
    overlap = {
        str(row["species"]): row
        for row in overlap_payload.get("species", [])
    }

    effort = {}
    with args.effort_csv.open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            if int(row["other_pilot_record_threshold"]) == effort_threshold:
                effort[str(row["species"])] = row

    ledger_payload = json.loads(
        args.occurrence_ledgers_json.read_text(encoding="utf-8")
    )
    if ledger_payload.get("schema") != LEDGER_MATRIX_SCHEMA:
        raise RuntimeError("unexpected occurrence-ledger matrix schema")
    transport = {}
    for ledger in ledger_payload.get("ledgers", []):
        for row in ledger.get("species", []):
            transport[str(row["species"])] = str(row["status"])

    rows = []
    qualified = []
    for name in species:
        reasons = []
        descriptor = descriptors.get(name)
        overlap_row = overlap.get(name)
        effort_row = effort.get(name)
        transport_status = transport.get(name)

        transport_ok = transport_status == "COMPLETE"
        if not transport_ok:
            reasons.append("OCCURRENCE_TRANSPORT_NOT_COMPLETE")

        if descriptor is None:
            host_families = None
            resolved_hosts = None
            host_ok = False
            reasons.append("MISSING_RESOURCE_DESCRIPTOR")
        else:
            host_families = float(descriptor["host_family_count"])
            resolved_hosts = int(descriptor["resolved_host_species"])
            host_ok = (
                host_families > 0
                and resolved_hosts >= host_families
            )
            if not host_ok:
                reasons.append("HOST_TAXONOMY_LOWER_BOUND_FAIL")

        if overlap_row is None:
            observed_units = 0
            outside = 0
            outside_fraction = None
            envelope_ok = False
            reasons.append("MISSING_CONTEMPORARY_OVERLAP")
        else:
            observed_units = int(overlap_row["butterfly_observed_units"])
            outside = int(overlap_row["observed_outside_contemporary"])
            outside_fraction = (
                None
                if observed_units == 0
                else outside / observed_units
            )
            envelope_ok = (
                outside_fraction is not None
                and outside_fraction <= maximum_outside
            )
            if not envelope_ok:
                reasons.append("CONTEMPORARY_ENVELOPE_OVERLAP_FAIL")

        if effort_row is None:
            supported = observed_supported = unobserved_supported = 0
            effort_ok = False
            reasons.append("MISSING_EFFORT_DIAGNOSTIC")
        else:
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

        passed = not reasons
        if passed:
            qualified.append(name)
        rows.append(
            {
                "species": name,
                "passed_preclimate_quality_gate": passed,
                "occurrence_transport_status": transport_status,
                "occurrence_transport_pass": transport_ok,
                "host_family_count": host_families,
                "resolved_host_species": resolved_hosts,
                "host_taxonomy_lower_bound_pass": host_ok,
                "butterfly_observed_units": observed_units,
                "outside_contemporary_units": outside,
                "outside_contemporary_fraction": outside_fraction,
                "contemporary_envelope_overlap_pass": envelope_ok,
                "effort_threshold_other_panel_records": effort_threshold,
                "effort_supported_host_units": supported,
                "effort_supported_observed_units": observed_supported,
                "effort_supported_unobserved_units": unobserved_supported,
                "sampling_effort_identifiability_pass": effort_ok,
                "failure_reasons": reasons,
            }
        )

    payload = {
        "schema": "ttf_butterfly_climate_release_preclimate_gate_v0.1",
        "status": (
            "PASS_TO_INDEPENDENT_CLIMATE_CROSSFIT"
            if len(qualified) >= minimum_species
            else str(q["if_below"])
        ),
        "panel_species_count": len(species),
        "panel_species_sha256": panel["species_sha256"],
        "qualified_species": qualified,
        "qualified_species_count": len(qualified),
        "minimum_species_required": minimum_species,
        "thresholds": thresholds,
        "species": rows,
        "climate_result_used": False,
        "pilot_climate_result_used_for_species_selection": False,
        "protocol_schema": protocol_schema,
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
