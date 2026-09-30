import json
import re
import runpy
from pathlib import Path

import pytest


ROOT=Path(__file__).resolve().parents[1]


def word_count(text: str) -> int:
    return len(text.split())


def test_ttf_q_main_document_meets_initial_mee_structure():
    text=(ROOT/"manuscript/ttf_q_methods_blinded_v0.1.md").read_text()
    match=re.search(
        r"## Abstract\n([\s\S]*?)\n## Data and code for peer review",
        text,
    )
    assert match is not None
    abstract=match.group(1)
    assert word_count(abstract) <= 350
    for number in ("1.","2.","3.","4."):
        assert number in abstract
    assert word_count(text) < 8000
    assert "## Figure legends" in text
    for number in range(1,6):
        assert f"Fig. {number}" in text
    assert "Table 1." in text
    assert "AI-assisted development disclosure" in text


def test_title_page_is_separate_and_contains_required_placeholders():
    main=(ROOT/"manuscript/ttf_q_methods_blinded_v0.1.md").read_text()
    title=(ROOT/"manuscript/ttf_q_title_page_template.md").read_text()
    assert "[AUTHOR 1 FULL NAME]" in title
    assert "[AFFILIATION 1]" in title
    assert "[EMAIL]" in title
    assert "Running headline" in title
    assert "Author contributions" in title
    assert "Conflict of interest" in title
    assert "[AUTHOR 1 FULL NAME]" not in main


def test_peer_review_bundle_is_fail_closed_without_author_license_choice(tmp_path):
    script=ROOT/"scripts/build_ttf_q_peer_review_bundle.py"
    namespace=runpy.run_path(str(script))
    assert namespace["ROOT"] == ROOT
    if (ROOT/"LICENSE").exists():
        pytest.skip("Author has now chosen a license; fail-closed no-license state no longer applies")
    import sys
    old_argv=sys.argv
    sys.argv=[str(script),"--output",str(tmp_path/"review.zip")]
    try:
        with pytest.raises(RuntimeError,match="LICENSE is required"):
            namespace["main"]()
    finally:
        sys.argv=old_argv


def test_submission_audit_does_not_claim_license_or_author_fields_complete():
    audit=json.loads(
        (ROOT/"manuscript/ttf_q_mee_submission_audit_v0.2.json").read_text()
    )
    assert audit["scientific_core"]["known_truth_information_validation"]=="COMPLETE"
    assert audit["scientific_core"]["genetic_response_opened"] is False
    assert audit["peer_review_code"]["license"]=="BLOCKING_AUTHOR_CHOICE"
    assert audit["title_page"]["author_names_affiliations_correspondence"]=="PENDING_AUTHOR_INPUT"



def test_anonymous_bundle_scanner_rejects_generic_identifying_metadata():
    namespace=runpy.run_path(
        str(ROOT/"scripts/build_ttf_q_peer_review_bundle.py")
    )
    scan=namespace["anonymity_hits"]
    assert scan("https://github.com/example/project")
    assert scan("contact: author@example.org")
    assert scan("https://orcid.org/0000-0000-0000-0000")
    assert not scan("https://doi.org/10.1111/example")


def test_new_prospective_stop_bundle_inputs_are_anonymous():
    namespace=runpy.run_path(
        str(ROOT/"scripts/build_ttf_q_peer_review_bundle.py")
    )
    scan=namespace["anonymity_hits"]
    paths=[
        "benchmarks/frozen/ttf_q_palearctic_prospective_stop_case_v0.1.json",
        "benchmarks/frozen/genetic_palearctic_lgm_v03_qualification_result_receipt_v0.1.json",
        "benchmarks/frozen/genetic_palearctic_lgm_v03_qualification_authority_correction_v0.1.json",
        "docs/supporting/genetic_palearctic_lgm_ttf_q_execution_v0.4.json",
    ]
    for rel in paths:
        body=(ROOT/rel).read_text()
        assert scan(body)==[], rel
