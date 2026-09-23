from pathlib import Path


QUALIFICATION = Path(".github/workflows/relational-environment-qualification-v03.yml")
POST = Path(".github/workflows/relational-environment-postqualification.yml")
HISTORICAL = Path(".github/workflows/relational-historical-prequalification.yml")


def test_relation_stage_not_evaluable_freezes_transition_receipt():
    text = QUALIFICATION.read_text()
    assert "Freeze relation-stage NOT_EVALUABLE_B transition" in text
    assert '"B_state":"NOT_EVALUABLE_B"' in text
    assert '"C_open_authorized":True' in text
    assert '"C_entry_state":"PROCEED_TO_C_AFTER_B_NOT_EVALUABLE"' in text
    assert "benchmarks/frozen/relational_environment_program_transition_v0.3.json" in text
    assert "git push origin HEAD:relational-ttf-v01" in text


def test_postqualification_persists_transition_to_repository():
    text = POST.read_text()
    assert "contents: write" in text
    assert "Commit frozen B-to-C transition receipt" in text
    assert "cp results/final/relational_environment_program_transition_v0.3.json" in text
    assert "benchmarks/frozen/relational_environment_program_transition_v0.3.json" in text
    assert "git push origin HEAD:relational-ttf-v01" in text


def test_study_c_reads_only_frozen_repository_transition():
    text = HISTORICAL.read_text()
    assert 'transition_rule["study_B_transition_receipt"]' in text
    assert "Study-C opening forbidden" in text
    assert "C_open_authorized" in text
