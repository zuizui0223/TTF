import runpy
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts/render_ttf_q_submission_files.py"
NS = runpy.run_path(str(SCRIPT))


def test_submission_renderer_reports_current_author_blockers():
    blockers = NS["readiness_blockers"]()
    if not (ROOT / "LICENSE").exists():
        assert "missing_root_LICENSE" in blockers
    title = ROOT / "manuscript/ttf_q_title_page_template.md"
    if title.exists() and any(
        token in title.read_text()
        for token in NS["TITLE_PLACEHOLDERS"]
    ):
        assert "unresolved_title_page_placeholders" in blockers


def test_renderer_requires_continuous_line_numbering_and_page_field_by_contract():
    source = SCRIPT.read_text()
    assert 'line_numbers.set(qn("w:countBy"), "1")' in source
    assert 'line_numbers.set(qn("w:restart"), "continuous")' in source
    assert 'instr.text = " PAGE "' in source
    assert "assert_submission_docx(main_docx, require_line_numbers=True)" in source
