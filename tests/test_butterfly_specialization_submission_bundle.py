from __future__ import annotations

import importlib.util
import re
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts" / "render_butterfly_specialization_blinded_manuscript.py"
SPEC = importlib.util.spec_from_file_location("render_blinded_butterfly_manuscript", SCRIPT)
assert SPEC is not None and SPEC.loader is not None
mod = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(mod)


def test_generated_blinded_manuscript_matches_renderer():
    source = (
        ROOT / "manuscript" / "butterfly_specialization_ecology_v0.1.md"
    ).read_text(encoding="utf-8")
    generated = (
        ROOT
        / "manuscript"
        / "generated"
        / "butterfly_specialization_ecology_blinded_v0.1.md"
    ).read_text(encoding="utf-8")
    assert generated == mod.render_blinded(source)


def test_blinded_manuscript_has_geb_required_front_matter():
    text = (
        ROOT
        / "manuscript"
        / "generated"
        / "butterfly_specialization_ecology_blinded_v0.1.md"
    ).read_text(encoding="utf-8")

    running = re.search(r"^\*\*Running title:\*\*\s*(.+)$", text, re.MULTILINE)
    assert running is not None
    assert len(running.group(1).strip()) < 40

    abstract_start = text.index("## Abstract")
    intro_start = text.index("## 1. Introduction")
    abstract = text[abstract_start:intro_start]
    for heading in (
        "Aim",
        "Location",
        "Time period",
        "Major taxa studied",
        "Methods",
        "Results",
        "Main conclusions",
    ):
        assert f"**{heading}:**" in abstract

    keywords = re.search(
        r"^\*\*Keywords:\*\*\s*(.+)$",
        abstract,
        re.MULTILINE,
    )
    assert keywords is not None
    items = [item.strip() for item in keywords.group(1).split(";") if item.strip()]
    if len(items) == 1:
        items = [item.strip() for item in keywords.group(1).split(",") if item.strip()]
    assert 6 <= len(items) <= 10
    assert items == sorted(items, key=str.casefold)


def test_blinded_manuscript_removes_internal_identity_tokens():
    text = (
        ROOT
        / "manuscript"
        / "generated"
        / "butterfly_specialization_ecology_blinded_v0.1.md"
    ).read_text(encoding="utf-8").lower()
    for token in (
        "ttf repository",
        ".github/workflows/",
        "repository provenance",
        "zuizui0223",
        "original ttf program",
    ):
        assert token not in text


def test_separate_title_page_template_contains_submission_metadata_slots():
    text = (
        ROOT
        / "manuscript"
        / "butterfly_specialization_geb_title_page_template_v0.1.md"
    ).read_text(encoding="utf-8")
    for required in (
        "Authors",
        "Affiliations",
        "Corresponding author",
        "Author contributions",
        "Acknowledgements",
        "Funding",
        "Conflict of interest",
        "Data and code availability",
        "Double-anonymous review note",
    ):
        assert required in text


def test_submission_checklist_points_to_double_anonymous_bundle():
    import json

    payload = json.loads(
        (
            ROOT
            / "manuscript"
            / "butterfly_specialization_geb_submission_checklist_v0.1.json"
        ).read_text(encoding="utf-8")
    )
    assert payload["journal"] == "Global Ecology and Biogeography"
    assert payload["author_guideline_snapshot"]["double_anonymous_review"] is True
    assert payload["author_guideline_snapshot"]["separate_title_page_required"] is True
    assert payload["manuscript"]["keyword_count"] == 7
    assert payload["figures"]["vector_pdf_available_for_all_five_figures"] is True
