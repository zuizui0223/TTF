from __future__ import annotations

import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
FAILURE = ROOT / "benchmarks" / "frozen" / "relational_historical_occurrence_terminal_failure_v0.1.json"
FINAL = ROOT / "benchmarks" / "frozen" / "relational_program_final_state_v0.3.json"
B_TRANSITION = ROOT / "benchmarks" / "frozen" / "relational_environment_program_transition_v0.3.json"


def _load(path: Path) -> dict:
    return json.loads(path.read_text())


def test_study_c_terminal_failure_preserves_frozen_contract() -> None:
    receipt = _load(FAILURE)

    assert receipt["schema"] == "ttf_relational_historical_occurrence_terminal_failure_v0.1"
    assert receipt["status"] == "IRREVERSIBLE_TECHNICAL_NOT_EVALUABLE_C_AFTER_FROZEN_RETRY_IMPLEMENTATION_FAILURE"

    run = receipt["workflow_run"]
    assert run["id"] == 35941577015
    assert run["head_sha"] == "ffacbd51d58689a4b18f7a2cb920f5a1385a74ab"
    assert run["conclusion"] == "failure"

    initial = receipt["initial_acquisition"]
    assert initial["required_batches"] == initial["completed_batches"] == 250
    assert initial["acquire_job_failures"] == 0
    assert initial["candidate_species"] == 1000
    assert sum(initial["status_counts_from_frozen_retry_plan"].values()) == 1000
    assert initial["status_counts_from_frozen_retry_plan"]["REQUEST_ERROR"] == 95

    plan = receipt["frozen_retry_plan"]
    assert plan["artifact_id"] == 10878157732
    assert plan["artifact_digest"] == "sha256:b4291a54514837f944f1cb9dbf14d320fd52293d24aa33b26d36ac063b7c7103"
    assert plan["request_error_count"] == 95
    assert plan["retry_batch_count"] == 24
    assert plan["maximum_retry_rounds"] == 1
    assert plan["scientific_query_change"] is False

    failure = receipt["observed_execution_failure"]
    assert failure["retry_jobs_required"] == failure["retry_jobs_failed"] == 24
    assert len(failure["failed_job_ids_by_batch"]) == 24
    assert failure["failure_occurs_before_batch_specific_retry_logic"] is True
    assert failure["retry_artifacts_uploaded"] == 0
    assert failure["finalize_job"] == "skipped"
    assert failure["final_artifact_present"] is False
    assert failure["frozen_binding_present"] is False
    assert failure["local_handoff_present"] is False

    unavailable = receipt["unavailable_final_fields"]
    for key in (
        "occurrence_status",
        "species_passing_ge_30",
        "final_request_error_count",
        "final_artifact_id",
        "final_artifact_digest",
    ):
        assert unavailable[key] is None

    decision = receipt["decision"]
    assert decision["C_state"] == "NOT_EVALUABLE_C"
    assert decision["biological_null"] is False
    assert decision["ecological_inadmissibility"] is False
    assert decision["replacement_run_authorized"] is False
    assert decision["rerun_authorized"] is False
    assert decision["threshold_change_authorized"] is False
    assert decision["downstream_Study_C_opening_authorized"] is False
    assert all(value is False for value in receipt["response_firewall"].values())


def test_terminal_failure_closes_finite_v03_without_biological_null() -> None:
    b = _load(B_TRANSITION)
    receipt = _load(FAILURE)
    final = _load(FINAL)

    assert b["B_state"] == "NOT_EVALUABLE_B"
    assert receipt["decision"]["C_state"] == "NOT_EVALUABLE_C"

    assert final == {
        "schema": "ttf_relational_program_final_state_v0.3",
        "B_state": "NOT_EVALUABLE_B",
        "C_state": "NOT_EVALUABLE",
        "C_decision": None,
        "C_p_one_sided": None,
        "final_state": "CLOSED_NO_EVALUABLE_TEST",
        "no_additional_predictor_authorized": True,
    }
