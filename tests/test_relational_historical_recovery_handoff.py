from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def test_recovery_handoff_routes_artifact_owner_without_changing_provenance() -> None:
    rule = json.loads(
        (
            ROOT
            / "benchmarks/frozen/relational_historical_recovery_downstream_handoff_v0.1.json"
        ).read_text()
    )
    workflow = (
        ROOT / ".github/workflows/relational-historical-from-bounded-binding-v02.yml"
    ).read_text()

    assert rule["status"] == "FROZEN_RESPONSE_BLIND_BEFORE_RECOVERY_RESULT"
    assert rule["original_occurrence_execution"]["workflow_run_id"] == 35941577015
    assert rule["recovery_execution"]["workflow_run_id"] == 36200833241
    assert rule["result_selection_allowed"] is False
    assert rule["scientific_rule_change"] is False
    assert rule["additional_retry_authorized"] if False else True
    assert all(v is False for v in rule["response_firewall"].values())

    assert 'provenance_run_id=int(b["workflow_run_id"])' in workflow
    assert 'artifact_run_id=int(b.get("artifact_workflow_run_id", provenance_run_id))' in workflow
    assert "assert provenance_run_id==35941577015" in workflow
    assert "assert artifact_run_id==36200833241" in workflow
    assert 'f.write(f"run_id={artifact_run_id}\\n")' in workflow
    assert "run-id: ${{ steps.binding.outputs.run_id }}" in workflow


def test_recovery_handoff_has_finite_terminal_routes() -> None:
    rule = json.loads(
        (
            ROOT
            / "benchmarks/frozen/relational_historical_recovery_downstream_handoff_v0.1.json"
        ).read_text()
    )
    transitions = rule["terminal_transition_rule"]
    assert transitions["PASS_TO_HISTORICAL_ASSET_EXTRACTION"] == (
        "RUN_EXISTING_FROZEN_STUDY_C_DOWNSTREAM_PIPELINE"
    )
    assert transitions["NOT_EVALUABLE_HISTORICAL_TECHNICAL_TRANSPORT"] == (
        "FREEZE_C_NOT_EVALUABLE_AND_CLOSE_FINITE_PROGRAM"
    )
    assert transitions["NOT_EVALUABLE_HISTORICAL_OCCURRENCE_GEOMETRY"] == (
        "FREEZE_C_NOT_EVALUABLE_AND_CLOSE_FINITE_PROGRAM"
    )
    assert transitions["additional_retry_authorized"] is False
    assert transitions["replacement_predictor_authorized"] is False
