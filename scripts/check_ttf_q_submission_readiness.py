#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import re
import runpy
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
MAIN = ROOT / "manuscript/ttf_q_methods_v0.1.md"
BLINDED = ROOT / "manuscript/ttf_q_methods_blinded_v0.1.md"
TITLE = ROOT / "manuscript/ttf_q_title_page_template.md"
ENQUIRY = ROOT / "manuscript/ttf_q_presubmission_enquiry_draft.md"
LICENSE = ROOT / "LICENSE"
BUNDLE_BUILDER = ROOT / "scripts/build_ttf_q_peer_review_bundle.py"
RENDERER = ROOT / "scripts/render_ttf_q_submission_files.py"
BLIND_BUILDER = ROOT / "scripts/build_ttf_q_blinded_manuscript.py"
FIGURE_RENDERER = ROOT / "scripts/render_ttf_q_figures.py"
PROSPECTIVE_STOP = ROOT / "benchmarks/frozen/ttf_q_palearctic_prospective_stop_case_v0.1.json"
SUBMISSION_AUDIT = ROOT / "manuscript/ttf_q_mee_submission_audit_v0.2.json"

TITLE_PLACEHOLDERS = (
    "[AUTHOR",
    "[AFFILIATION",
    "[INSTITUTION",
    "[NAME]",
    "[POSTAL ADDRESS]",
    "[EMAIL]",
    "[COMPLETE FOR TITLE PAGE",
    "[CRediT",
    "[DECLARATION]",
)


def words(text: str) -> int:
    return len(text.split())


def _blinded_expected(source: str) -> str:
    namespace = runpy.run_path(str(BLIND_BUILDER))
    return namespace["blind"](source)


def check_state(root: Path = ROOT) -> dict[str, object]:
    # Root is kept as an argument for tests/future packaging, but repository
    # paths above are authoritative for the checked-in TTF-Q submission state.
    del root

    structural_failures: list[str] = []
    hard_blockers: list[str] = []
    author_actions: list[str] = []
    non_blocking: list[str] = []

    required = [
        MAIN,
        BLINDED,
        TITLE,
        ENQUIRY,
        BUNDLE_BUILDER,
        RENDERER,
        BLIND_BUILDER,
        ROOT / "docs/TTF_Q_THEORY_V01.md",
        ROOT / "benchmarks/frozen/ttf_q_known_truth_detectability_v0.2.json",
        ROOT / "benchmarks/frozen/ttf_q_bc_calibrated_detectability_result_v0.2.json",
        FIGURE_RENDERER,
        PROSPECTIVE_STOP,
        SUBMISSION_AUDIT,
    ]
    missing = [str(path.relative_to(ROOT)) for path in required if not path.is_file()]
    if missing:
        structural_failures.append("missing_required_files:" + ",".join(missing))

    main_text = MAIN.read_text() if MAIN.is_file() else ""
    blinded_text = BLINDED.read_text() if BLINDED.is_file() else ""

    abstract_match = re.search(
        r"## Abstract\n([\s\S]*?)\n## Data and code for peer review",
        blinded_text,
    )
    abstract_words = words(abstract_match.group(1)) if abstract_match else None
    full_words = words(blinded_text) if blinded_text else None

    if abstract_match is None:
        structural_failures.append("missing_numbered_abstract")
    else:
        for number in ("1.", "2.", "3.", "4."):
            if number not in abstract_match.group(1):
                structural_failures.append(f"abstract_missing_item_{number[0]}")
        if abstract_words is not None and abstract_words > 350:
            structural_failures.append("abstract_over_350_words")

    if full_words is not None and full_words >= 8000:
        structural_failures.append("main_document_over_8000_words")

    for number in range(1, 6):
        if f"Fig. {number}" not in blinded_text:
            structural_failures.append(f"missing_figure_reference_{number}")
    if "## Figure legends" not in blinded_text:
        structural_failures.append("missing_figure_legends")
    if "Table 1." not in blinded_text:
        structural_failures.append("missing_table_1")
    if "AI-assisted development disclosure" not in blinded_text:
        structural_failures.append("missing_AI_disclosure")

    if PROSPECTIVE_STOP.is_file():
        case=json.loads(PROSPECTIVE_STOP.read_text())
        if case.get("status")!="RESPONSE_BLIND_STOP_BEFORE_SUBPANEL_GENETIC_RESPONSE":
            structural_failures.append("prospective_stop_receipt_status_drift")
        firewall=case.get("biological_response_firewall",{})
        if any(bool(v) for v in firewall.values()):
            structural_failures.append("prospective_stop_genetic_firewall_open")
        det=case.get("authoritative_deterministic_ttf_q",{})
        if det.get("overall_pass") is not False:
            structural_failures.append("prospective_stop_decision_drift")

    if SUBMISSION_AUDIT.is_file():
        audit=json.loads(SUBMISSION_AUDIT.read_text())
        core=audit.get("scientific_core",{})
        if core.get("prospective_response_blind_stop_application")!="COMPLETE":
            structural_failures.append("submission_audit_missing_prospective_stop")
        if core.get("figures_from_frozen_scalar_handoff")!="COMPLETE_FIVE_FIGURES":
            structural_failures.append("submission_audit_figure_count_drift")

    if FIGURE_RENDERER.is_file():
        figure_text=FIGURE_RENDERER.read_text()
        if "figure5_prospective_stop.svg" not in figure_text:
            structural_failures.append("figure5_renderer_missing")

    if MAIN.is_file() and BLINDED.is_file() and BLIND_BUILDER.is_file():
        if blinded_text != _blinded_expected(main_text):
            structural_failures.append("blinded_document_not_deterministic_from_source")

    lower_blinded = blinded_text.lower()
    for token in ("github.com/zuizui0223", "zuizui0223/ttf"):
        if token in lower_blinded:
            structural_failures.append("public_repository_identity_in_blinded_document")

    if not LICENSE.is_file():
        hard_blockers.append("author_must_choose_open_source_license")

    title_text = TITLE.read_text() if TITLE.is_file() else ""
    unresolved = sorted(token for token in TITLE_PLACEHOLDERS if token in title_text)
    if unresolved:
        hard_blockers.append("author_must_complete_title_page_fields")

    if ENQUIRY.is_file() and "[CORRESPONDING AUTHOR]" in ENQUIRY.read_text():
        author_actions.append("send_or_decline_presubmission_scope_enquiry")

    if not (ROOT / "CITATION.cff").is_file():
        non_blocking.append("final_CITATION_cff_pending")
    non_blocking.append("permanent_archive_DOI_pending_until_public_release")

    scientific_ready = not structural_failures
    final_render_ready = scientific_ready and not hard_blockers
    if not scientific_ready:
        status = "SCIENTIFIC_OR_PACKAGING_STRUCTURE_ERROR"
    elif hard_blockers:
        status = "AUTHOR_INPUT_REQUIRED"
    else:
        status = "READY_FOR_FINAL_RENDER"

    return {
        "schema": "ttf_q_submission_readiness_check_v0.1",
        "status": status,
        "scientific_package_ready": scientific_ready,
        "final_render_ready": final_render_ready,
        "structural_failures": structural_failures,
        "hard_blockers": hard_blockers,
        "author_actions": author_actions,
        "non_blocking": non_blocking,
        "metrics": {
            "abstract_words": abstract_words,
            "blinded_document_words": full_words,
            "abstract_limit": 350,
            "main_document_limit": 8000,
            "expected_figure_count": 5,
        },
        "expected_next_step": (
            "Fix structural failures before submission work."
            if structural_failures
            else (
                "Resolve author license/title-page choices, then run final renderer."
                if hard_blockers
                else "Run final renderer and visually audit the generated submission files."
            )
        ),
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path)
    parser.add_argument(
        "--require-ready",
        action="store_true",
        help="Exit non-zero unless all author blockers are resolved.",
    )
    args = parser.parse_args()

    payload = check_state()
    rendered = json.dumps(payload, indent=2, sort_keys=True) + "\n"
    if args.output is not None:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered)
    print(rendered, end="")

    if payload["structural_failures"]:
        return 2
    if args.require_ready and payload["hard_blockers"]:
        return 3
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
