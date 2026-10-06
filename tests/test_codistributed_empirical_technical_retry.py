from pathlib import Path


ROOT=Path(__file__).resolve().parents[1]


def test_codistributed_empirical_technical_continuation_is_predeclared_and_guarded():
    rule=(ROOT/"benchmarks/frozen/genetic_codistributed_recurrence_empirical_technical_retry_rule_v0.1.json").read_text()
    workflow=(ROOT/".github/workflows/codistributed-recurrence-empirical-technical-continuation.yml").read_text()
    assert '"status": "FROZEN_BEFORE_FIRST_EMPIRICAL_RUN_START"' in rule
    assert '"workflow_run_id": 37424857802' in rule
    assert "workflow_dispatch" not in workflow
    assert 'FIRST_EMPIRICAL_RUN_ID: "37424857802"' in workflow
    assert 'FIRST_EMPIRICAL_JOB_ID: "112142102812"' in workflow
    assert 'assert empirical["conclusion"]=="skipped"' in workflow
    assert "fail-on-cache-miss: true" in workflow
    assert "cmp \"$AUTHORIZATION\" results/pre-opening/reproduced_authorization.json" in workflow
    assert "run_codistributed_recurrence_confirmatory_empirical.py" in workflow
