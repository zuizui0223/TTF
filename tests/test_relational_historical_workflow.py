from pathlib import Path


def test_historical_workflow_is_manual_and_transition_guarded():
    text = Path(".github/workflows/relational-historical-prequalification.yml").read_text()
    assert "workflow_dispatch:" in text
    assert "push:" not in text
    assert "relational_environment_program_transition_v0.3.json" not in text
    assert 'transition_rule["study_B_transition_receipt"]' in text
    assert 'C_open_authorized' in text
    assert 'DETECTED_B' not in text  # decision is read from frozen transition contract, not reimplemented here


def test_historical_workflow_stops_before_genetic_response():
    text = Path(".github/workflows/relational-historical-prequalification.yml").read_text()
    assert "run_relational_historical_climate_qualification.py" in text
    assert "PASS_TO_C_CONFIRMATORY_CHARACTER_MASK" in text
    assert "Study_C_genetic_response_opened" in text
    assert "run_relational" in text
    assert "empirical" not in text.lower()
