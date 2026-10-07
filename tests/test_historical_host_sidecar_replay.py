from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]

SIDECAR_RUN="35548477097"
INTERACTION_SHA="0a084fb5273e4780b03e015205c87e4545c69f0f8e517feb1a2a3879889508e9"
NATIVE_SHA="c731906315f7452f83ce302c1d36c75ad94afc24ce7e24d461360312244ac558"


def test_occurrence_workflow_replays_exact_sidecar():
    text=(ROOT/".github/workflows/historical-host-occurrence-gate.yml").read_text()
    assert f"run-id: {SIDECAR_RUN}" in text
    assert INTERACTION_SHA in text
    assert NATIVE_SHA in text
    assert "export_historical_host_native_ranges_from_sidecar.R" in text
    assert "_host_sidecar/native_extant_nondoubtful_wgsrpd3.csv" in text


def test_neotoma_workflow_replays_exact_sidecar():
    text=(ROOT/".github/workflows/historical-host-neotoma-validation.yml").read_text()
    assert f"run-id: {SIDECAR_RUN}" in text
    assert NATIVE_SHA in text
    assert "export_historical_host_names_from_sidecar.py" in text
    assert "_host_sidecar/native_extant_nondoubtful_wgsrpd3.csv" in text


def test_sidecar_repair_keeps_scientific_invariants():
    import json
    p=json.loads((ROOT/"benchmarks/frozen/genetic_historical_host_sidecar_replay_repair_v0.1.json").read_text())
    assert p["status"]=="FROZEN_AFTER_PREQUERY_ID_RECONSTRUCTION_FAILURE_BEFORE_ANY_AUTHORITATIVE_EXTERNAL_DATA_RESULT"
    assert p["scientific_invariants"]["insect_species"]==140
    assert p["scientific_invariants"]["host_ids_changed"] is False
    assert p["scientific_invariants"]["insect_membership_changed"] is False
    assert p["scientific_invariants"]["gbif_minimum_occurrences"]==30
    assert p["scientific_invariants"]["insect_host_completeness"]==0.8
    assert all(v is False for v in p["response_firewall"].values())


def test_replacement_run_authority_is_pre_result_and_exact():
    import json
    p=json.loads((ROOT/"benchmarks/frozen/genetic_historical_host_external_execution_authority_v0.2.json").read_text())
    assert p["status"]=="FROZEN_BEFORE_SIDECAR_REPLAY_REPLACEMENT_RESULTS"
    assert p["replacement_authoritative_runs"]["native_range_host_gbif"]["workflow_run_id"]==37644277020
    assert p["replacement_authoritative_runs"]["neotoma_lgm_validation"]["workflow_run_id"]==37644319101
    assert p["replacement_authoritative_runs"]["native_range_host_gbif"]["result_seen_at_freeze"] is False
    assert p["replacement_authoritative_runs"]["neotoma_lgm_validation"]["result_seen_at_freeze"] is False
