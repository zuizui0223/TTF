from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def test_historical_qualification_runner_has_survivor_stage():
    text = (ROOT / "scripts/run_relational_historical_climate_qualification.py").read_text()
    assert '--qualification-stage' in text
    assert 'choices=("development", "confirmatory_survivor")' in text
    assert '"PASS_TO_HISTORICAL_EMPIRICAL_IDENTITY_OPENING_PREPARATION"' in text
    assert '"NOT_EVALUABLE_C_CHARACTER_MASK_OR_SURVIVOR_GEOMETRY"' in text
    assert '"qualification_stage": args.qualification_stage' in text


def test_historical_workflow_is_gated_by_B_transition_and_closes_program():
    text = (ROOT / ".github/workflows/relational-historical-prequalification.yml").read_text()
    assert "Study-C opening forbidden" in text
    assert "--qualification-stage confirmatory_survivor" in text
    assert "freeze_relational_historical_empirical_authorization.py" in text
    assert "run_relational_historical_empirical.py" in text
    assert '"DETECTED_C_AFTER_B_NULL"' in text
    assert '"DETECTED_C_AFTER_B_NOT_EVALUABLE"' in text
    assert '"CLOSED_TWO_NULLS"' in text
    assert '"CLOSED_PARTIAL_NOT_EVALUABLE"' in text
    assert '"CLOSED_NO_EVALUABLE_TEST"' in text
