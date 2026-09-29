import json
import runpy
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts/check_ttf_q_submission_readiness.py"
NS = runpy.run_path(str(SCRIPT))


def test_current_submission_state_has_no_scientific_structure_failure():
    state = NS["check_state"]()
    assert state["scientific_package_ready"] is True
    assert state["structural_failures"] == []
    assert state["metrics"]["abstract_words"] <= 350
    assert state["metrics"]["blinded_document_words"] < 8000


def test_current_remaining_hard_blockers_are_author_choices_only():
    state = NS["check_state"]()
    expected = set()
    if not (ROOT / "LICENSE").exists():
        expected.add("author_must_choose_open_source_license")

    title = ROOT / "manuscript/ttf_q_title_page_template.md"
    placeholders = NS["TITLE_PLACEHOLDERS"]
    if title.exists() and any(token in title.read_text() for token in placeholders):
        expected.add("author_must_complete_title_page_fields")

    assert set(state["hard_blockers"]) == expected
    if expected:
        assert state["status"] == "AUTHOR_INPUT_REQUIRED"
        assert state["final_render_ready"] is False


def test_readiness_checker_identifies_presubmission_enquiry_as_author_action_not_scientific_blocker():
    state = NS["check_state"]()
    assert "send_or_decline_presubmission_scope_enquiry" in state["author_actions"]
    assert "send_or_decline_presubmission_scope_enquiry" not in state["hard_blockers"]


def test_readiness_payload_is_json_serializable():
    json.dumps(NS["check_state"]())
