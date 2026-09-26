from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def test_study_c_retry_terminal_inevitability_receipt_matches_frozen_finalizer() -> None:
    receipt = json.loads(
        (
            ROOT
            / "benchmarks/frozen/relational_historical_retry_terminal_inevitability_v0.1.json"
        ).read_text()
    )
    recovery = json.loads(
        (
            ROOT
            / "benchmarks/frozen/relational_historical_retry_import_recovery_v0.1.json"
        ).read_text()
    )
    finalizer = (
        ROOT / "scripts/finalize_relational_historical_occurrence_bounded.py"
    ).read_text()

    assert receipt["status"] == (
        "OBSERVED_IRREVERSIBLE_OCCURRENCE_CLASS_IF_FROZEN_FINALIZER_COMPLETES"
    )
    assert receipt["recovery_workflow"]["run_id"] == 36200833241
    assert receipt["recovery_workflow"]["retry_round"] == 1
    assert receipt["recovery_workflow"]["additional_retry_authorized"] is False
    assert recovery["recovery_execution"]["additional_retry_round_authorized"] is False

    statuses = []
    for batch in receipt["completed_retry_batches_observed"]:
        assert batch["job_conclusion"] == "success"
        statuses.extend(batch["species_statuses"].values())
    assert len(statuses) == receipt["observed_retry_species"] == 8
    assert statuses.count("REQUEST_ERROR") == receipt["observed_request_error_after_retry"] == 8

    assert 'unresolved=[name for name in names if species[name]["status"]=="REQUEST_ERROR"]' in finalizer
    assert '"NOT_EVALUABLE_HISTORICAL_TECHNICAL_TRANSPORT"' in finalizer
    assert receipt["frozen_finalizer_logic"][
        "classification_if_exact_finalizer_runs_with_these_artifacts"
    ] == "NOT_EVALUABLE_HISTORICAL_TECHNICAL_TRANSPORT"

    interpretation = receipt["interpretation"]
    assert interpretation[
        "remaining_retry_outcomes_can_no_longer_restore_PASS_TO_HISTORICAL_ASSET_EXTRACTION"
    ] is True
    assert interpretation["remaining_batches_must_still_complete_for_exact_recovery_audit"] is True
    assert interpretation["early_workflow_cancellation_authorized"] is False
    assert interpretation["canonical_final_state_may_be_written_before_recovery_result"] is False
    assert interpretation["biological_null"] is False
    assert receipt["relation_result_seen"] is False
    assert receipt["genetic_response_used"] is False
    assert all(v is False for v in receipt["response_firewall"].values())
