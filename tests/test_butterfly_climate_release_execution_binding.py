from __future__ import annotations

import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def test_independent_climate_release_execution_binding_matches_frozen_contract():
    binding = json.loads(
        (
            ROOT
            / "benchmarks/exploratory/butterfly_climate_release_execution_binding_v0.1.json"
        ).read_text()
    )
    protocol = json.loads(
        (
            ROOT
            / "docs/exploratory/butterfly_climate_release_independent_test_v0.1.json"
        ).read_text()
    )
    recovery = json.loads(
        (
            ROOT
            / "docs/exploratory/butterfly_climate_release_transport_recovery_v0.1.json"
        ).read_text()
    )
    panel = json.loads(
        (
            ROOT
            / "benchmarks/exploratory/butterfly_climate_release_independent_panel_v0.1.json"
        ).read_text()
    )

    assert binding["status"] == "SUPERSEDED_BEFORE_PRECLIMATE_OR_CLIMATE_RESULT"
    assert binding["execution"]["run_id"] == 36239130041
    assert binding["execution"]["head_sha"] == (
        "4c2fd1093acfb72c9b9c7b681c46861a6a29aae3"
    )

    assert panel["panel"]["species_count"] == binding["independent_panel"]["species_count"] == 32
    assert panel["panel"]["species_sha256"] == binding["independent_panel"]["species_sha256"]
    assert panel["panel"]["species_sha256"] == (
        "bbaaa28cf361fb02d558988bdd1a3f948a243b79c8d3bb1709e908edb1f0f204"
    )

    primary = protocol["primary_test"]
    bound = binding["scientific_contract"]
    assert primary["alternative"] == bound["primary_alternative"] == "negative"
    assert int(primary["permutation"]["iterations"]) == bound["permutations"] == 9999
    assert float(primary["alpha"]) == bound["alpha"] == 0.05
    assert int(
        protocol["preclimate_quality_gate"][
            "minimum_species_passing_preclimate_gate"
        ]
    ) == bound["minimum_preclimate_qualified_species"] == 12
    assert int(
        protocol["climate_crossfit"]["minimum_climate_informative_species"]
    ) == bound["minimum_climate_informative_species"] == 12

    technical = recovery["technical_change"]
    transport = binding["transport_contract"]
    assert int(technical["original_window_size_records"]) == (
        transport["ordinal_window_size_records"]
    ) == 300
    assert int(technical["recovery_transport_chunk_size_records"]) == (
        transport["recovery_chunk_size_records"]
    ) == 10
    assert int(technical["request_timeout_cap_seconds"]) == (
        transport["request_timeout_seconds"]
    ) == 60
    assert int(technical["species_total_deadline_seconds_per_pass"]) == (
        transport["species_deadline_seconds_per_pass"]
    ) == 600
    assert int(technical["maximum_passes_per_job"]) == (
        transport["maximum_passes_per_job"]
    ) == 2

    assert all(
        value is False
        for value in binding["response_firewall_at_binding"].values()
    )


def test_v01_execution_was_superseded_without_opening_ecological_response():
    binding = json.loads(
        (
            ROOT
            / "benchmarks/exploratory/butterfly_climate_release_execution_binding_v0.1.json"
        ).read_text()
    )
    supersession = binding["supersession"]
    assert supersession["independent_host_overlap_opened_before_supersession"] is False
    assert supersession["preclimate_quality_gate_opened_before_supersession"] is False
    assert supersession["independent_climate_values_opened_before_supersession"] is False
    assert supersession["independent_climate_filtering_scores_opened_before_supersession"] is False
    assert supersession["independent_primary_result_opened_before_supersession"] is False
    assert supersession["transport_only_state_may_be_reused"] is True
    assert supersession["scientific_result_from_run_36239130041_must_not_be_used"] is True


def test_v021_execution_binding_matches_row_order_invariant_protocol():
    binding = json.loads(
        (
            ROOT
            / "benchmarks/exploratory/butterfly_climate_release_execution_binding_v0.2.1.json"
        ).read_text()
    )
    protocol = json.loads(
        (
            ROOT
            / "docs/exploratory/butterfly_climate_release_independent_test_v0.2.1.json"
        ).read_text()
    )
    panel = json.loads(
        (
            ROOT
            / "benchmarks/exploratory/butterfly_climate_release_independent_panel_v0.1.json"
        ).read_text()
    )

    assert binding["status"] == "BOUND_BEFORE_PRECLIMATE_OR_CLIMATE_RESULT"
    assert binding["execution"]["run_id"] == 36241187187
    assert binding["execution"]["head_sha"] == (
        "0d4b19a67652dd233c0632eff9deadfb865bc7c3"
    )
    assert binding["scientific_contract"]["protocol_schema"] == protocol["schema"]
    assert protocol["schema"] == (
        "ttf_butterfly_climate_release_independent_test_v0.2.1"
    )
    assert protocol["primary_test"]["permutation"]["iterations"] == 9999
    assert protocol["primary_test"]["alternative"] == "negative"
    assert protocol["primary_test"]["alpha"] == 0.05
    assert "species-identity" in protocol["primary_test"]["permutation"]["method"]
    assert (
        panel["panel"]["species_sha256"]
        == binding["independent_panel"]["species_sha256"]
    )
    assert all(
        value is False
        for value in binding["response_firewall_at_binding"].values()
    )
