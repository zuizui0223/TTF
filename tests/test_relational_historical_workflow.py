from pathlib import Path


WORKFLOW = Path(".github/workflows/relational-historical-prequalification.yml")


def test_historical_workflow_is_transition_triggered_and_guarded():
    text = WORKFLOW.read_text()
    assert "workflow_dispatch:" in text
    assert "push:" in text
    assert 'benchmarks/frozen/relational_environment_program_transition_v0.3.json' in text
    assert 'transition_rule["study_B_transition_receipt"]' in text
    assert "Resolve frozen B-to-C transition" in text
    assert "C_open_authorized" in text
    assert "Study-C opening forbidden" in text
    assert "needs.authorize.outputs.authorized == 'true'" in text
    assert 'assert state=="DETECTED_B"' in text
    assert 'assert receipt.get("C_open_authorized") is False' in text


def test_historical_workflow_contains_frozen_full_pipeline_after_guard():
    text = WORKFLOW.read_text()
    assert "freeze_relational_historical_candidate_census.py" in text
    assert "acquire_relational_historical_occurrences.py" in text
    assert "build_relational_historical_relation.py" in text
    assert "attach_relational_historical_opportunity.py" in text
    assert "run_relational_historical_climate_qualification.py" in text
    assert "freeze_relational_historical_character_mask.py" in text
    assert "--qualification-stage confirmatory_survivor" in text
    assert "freeze_relational_historical_empirical_authorization.py" in text
    assert "run_relational_historical_empirical.py" in text


def test_historical_workflow_closes_finite_program():
    text = WORKFLOW.read_text()
    assert "relational_program_final_state_v0.3.json" in text
    assert '"DETECTED_C_AFTER_B_NULL"' in text
    assert '"DETECTED_C_AFTER_B_NOT_EVALUABLE"' in text
    assert '"CLOSED_TWO_NULLS"' in text
    assert '"CLOSED_PARTIAL_NOT_EVALUABLE"' in text
    assert '"CLOSED_NO_EVALUABLE_TEST"' in text
    assert '"no_additional_predictor_authorized":True' in text


def test_historical_workflow_uses_canonical_phylogatr_marker_contract():
    text = WORKFLOW.read_text()
    assert "--protocol docs/supporting/genetic_phylogatr_confirmatory_protocol_v0.1.json" in text
    assert "--protocol docs/supporting/relational_genetic_protocol_v0.1.json" not in text
