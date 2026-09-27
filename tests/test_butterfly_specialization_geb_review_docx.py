from __future__ import annotations

import importlib.util
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts" / "render_butterfly_specialization_geb_review_docx.py"
SPEC = importlib.util.spec_from_file_location(
    "render_butterfly_specialization_geb_review_docx",
    SCRIPT,
)
assert SPEC is not None and SPEC.loader is not None
mod = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(mod)


def test_blinded_markdown_uses_geb_section_order():
    text = (
        ROOT
        / "manuscript"
        / "generated"
        / "butterfly_specialization_ecology_blinded_v0.1.md"
    ).read_text(encoding="utf-8")
    sections = mod.split_submission_sections(text)
    assert "## References" in sections["main_through_data"]
    assert "## Data and Code Availability" in sections["main_through_data"]
    assert sections["figure_legends"].startswith("## Figure legends")
    assert text.index("## References") < text.index("## Data and Code Availability")
    assert text.index("## Data and Code Availability") < text.index("## Figure legends")


def test_blinded_keywords_are_comma_separated_alphabetized_and_in_range():
    text = (
        ROOT
        / "manuscript"
        / "generated"
        / "butterfly_specialization_ecology_blinded_v0.1.md"
    ).read_text(encoding="utf-8")
    items = mod.keyword_items(text)
    assert 6 <= len(items) <= 10
    assert items == sorted(items, key=str.casefold)
    line = next(
        line for line in text.splitlines() if line.startswith("**Keywords:**")
    )
    assert ";" not in line
    assert "," in line


def test_exactly_five_figure_legends_are_available_for_docx_embedding():
    text = (
        ROOT
        / "manuscript"
        / "generated"
        / "butterfly_specialization_ecology_blinded_v0.1.md"
    ).read_text(encoding="utf-8")
    legends = mod.extract_figure_legends(
        text[text.index("## Figure legends") :]
    )
    assert sorted(legends) == [1, 2, 3, 4, 5]
    for number, legend in legends.items():
        assert legend.startswith(f"**Figure {number}.")
        assert len(legend) > 80
