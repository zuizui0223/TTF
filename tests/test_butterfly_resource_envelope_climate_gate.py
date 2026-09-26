from __future__ import annotations

import csv
import json
import subprocess
import sys
from pathlib import Path


SCRIPT = Path("scripts/apply_butterfly_resource_envelope_climate_gate.py")


def _write_csv(path: Path, fieldnames: list[str], rows: list[dict[str, object]]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def test_preclimate_gate_requires_complete_transport_and_quality(tmp_path: Path):
    gate = {
        "schema": "ttf_butterfly_resource_envelope_climate_pilot_gate_v0.1",
        "species_quality_gate": {
            "sampling_effort_identifiability": {
                "other_pilot_record_threshold": 10,
                "minimum_effort_supported_host_units": 5,
                "minimum_butterfly_observed_effort_supported_host_units": 2,
                "minimum_butterfly_unobserved_effort_supported_host_units": 2,
            },
            "contemporary_envelope_overlap": {
                "maximum_outside_contemporary_fraction": 0.10
            },
            "climate_training_floor": {
                "minimum_valid_training_occurrence_records": 30
            },
        },
        "pilot_level_gate": {
            "minimum_species_passing_quality_gate": 1,
            "if_below": "STOP_CLIMATE_PILOT_AS_NOT_IDENTIFIABLE_AT_CURRENT_RESOLUTION",
        },
    }
    gate_path = tmp_path / "gate.json"
    gate_path.write_text(json.dumps(gate), encoding="utf-8")

    descriptors = tmp_path / "descriptors.csv"
    _write_csv(
        descriptors,
        ["species", "host_family_count", "resolved_host_species"],
        [
            {"species": "Pass species", "host_family_count": 2, "resolved_host_species": 4},
            {"species": "Partial species", "host_family_count": 1, "resolved_host_species": 3},
            {"species": "Bad envelope", "host_family_count": 1, "resolved_host_species": 2},
        ],
    )

    overlap = {
        "schema": "ttf_butterfly_resource_envelope_native_contemporary_occurrence_overlap_v0.1",
        "species": [
            {
                "species": "Pass species",
                "butterfly_observed_units": 10,
                "observed_outside_contemporary": 1,
            },
            {
                "species": "Partial species",
                "butterfly_observed_units": 10,
                "observed_outside_contemporary": 0,
            },
            {
                "species": "Bad envelope",
                "butterfly_observed_units": 10,
                "observed_outside_contemporary": 3,
            },
        ],
    }
    overlap_path = tmp_path / "overlap.json"
    overlap_path.write_text(json.dumps(overlap), encoding="utf-8")

    effort = tmp_path / "effort.csv"
    _write_csv(
        effort,
        [
            "species",
            "other_pilot_record_threshold",
            "effort_supported_host_units",
            "occupied_effort_supported_host_units",
            "unoccupied_effort_supported_host_units",
        ],
        [
            {
                "species": name,
                "other_pilot_record_threshold": 10,
                "effort_supported_host_units": 8,
                "occupied_effort_supported_host_units": 3,
                "unoccupied_effort_supported_host_units": 5,
            }
            for name in ("Pass species", "Partial species", "Bad envelope")
        ],
    )

    ledgers = {
        "schema": "ttf_butterfly_resource_envelope_occurrence_matrix_v0.1",
        "ledgers": [
            {"species": [{"species": "Pass species", "status": "COMPLETE"}]},
            {"species": [{"species": "Partial species", "status": "PARTIAL_TRANSPORT"}]},
            {"species": [{"species": "Bad envelope", "status": "COMPLETE"}]},
        ],
    }
    ledgers_path = tmp_path / "ledgers.json"
    ledgers_path.write_text(json.dumps(ledgers), encoding="utf-8")

    output = tmp_path / "result.json"
    subprocess.run(
        [
            sys.executable,
            str(SCRIPT),
            "--gate-json",
            str(gate_path),
            "--descriptors-csv",
            str(descriptors),
            "--overlap-json",
            str(overlap_path),
            "--effort-csv",
            str(effort),
            "--occurrence-ledgers-json",
            str(ledgers_path),
            "--output-json",
            str(output),
        ],
        check=True,
        capture_output=True,
        text=True,
    )
    result = json.loads(output.read_text(encoding="utf-8"))
    assert result["status"] == "PASS_TO_EXPLORATORY_CLIMATE_PILOT"
    assert result["qualified_species"] == ["Pass species"]

    rows = {row["species"]: row for row in result["species"]}
    assert rows["Pass species"]["occurrence_transport_pass"] is True
    assert rows["Partial species"]["occurrence_transport_pass"] is False
    assert "OCCURRENCE_TRANSPORT_NOT_COMPLETE" in rows["Partial species"]["failure_reasons"]
    assert "CONTEMPORARY_ENVELOPE_OVERLAP_FAIL" in rows["Bad envelope"]["failure_reasons"]
