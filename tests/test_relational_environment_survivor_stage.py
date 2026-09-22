from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def test_environment_qualification_runner_has_explicit_survivor_stage():
    text = (ROOT / "scripts/run_relational_environment_qualification.py").read_text()
    assert '--qualification-stage' in text
    assert 'choices=("development", "confirmatory_survivor")' in text
    assert '"PASS_TO_EMPIRICAL_IDENTITY_OPENING_PREPARATION"' in text
    assert '"NOT_EVALUABLE_B_CHARACTER_MASK_OR_SURVIVOR_GEOMETRY"' in text
    assert '"qualification_stage": args.qualification_stage' in text


def test_postqualification_workflow_requests_survivor_stage():
    text = (ROOT / ".github/workflows/relational-environment-postqualification.yml").read_text()
    assert "--qualification-stage confirmatory_survivor" in text
