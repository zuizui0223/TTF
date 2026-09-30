#!/usr/bin/env python3
from __future__ import annotations

import argparse
import hashlib
import json
import re
import shutil
import tempfile
import zipfile
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]

FILES = [
    "LICENSE",
    "pyproject.toml",
    "manuscript/ttf_q_methods_blinded_v0.1.md",
    "manuscript/generated/ttf_q_scalar_results_v0.1.json",
    "docs/TTF_Q_V01.md",
    "docs/TTF_Q_THEORY_V01.md",
    "docs/supporting/ttf_q_v0.1.json",
    "docs/supporting/ttf_q_known_truth_benchmark_v0.1.json",
    "docs/supporting/ttf_q_known_truth_detectability_v0.2.json",
    "benchmarks/frozen/ttf_q_c1_result_receipt_v0.1.json",
    "benchmarks/frozen/ttf_q_c2_result_receipt_v0.1.json",
    "benchmarks/frozen/ttf_q_b1_c1_information_ladder_v0.1.json",
    "benchmarks/frozen/ttf_q_bc_calibrated_detectability_result_v0.2.json",
    "benchmarks/frozen/ttf_q_known_truth_decomposition_v0.1.json",
    "benchmarks/frozen/ttf_q_known_truth_detectability_v0.2.json",
    "benchmarks/frozen/ttf_q_palearctic_prospective_stop_case_v0.1.json",
    "benchmarks/frozen/genetic_palearctic_lgm_v03_qualification_result_receipt_v0.1.json",
    "benchmarks/frozen/genetic_palearctic_lgm_v03_qualification_authority_correction_v0.1.json",
    "docs/supporting/genetic_palearctic_lgm_ttf_q_execution_v0.4.json",
    "src/ttf/relational_dyadic.py",
    "src/ttf/relational_qualification.py",
    "src/ttf/relational_benchmark.py",
    "src/ttf/ecological_information_ladder.py",
    "src/ttf/external_geographic_opportunity.py",
    "src/ttf/historical_geographic_information.py",
    "src/ttf/historical_relation_information.py",
    "src/ttf/climate_occurrence_filter.py",
    "src/ttf/relational_environment.py",
    "src/ttf/relational_historical.py",
    "scripts/run_ttf_q_characterization.py",
    "scripts/run_ttf_q_known_truth_benchmark.py",
    "scripts/run_ttf_q_known_truth_detectability.py",
    "scripts/aggregate_ttf_q_known_truth_detectability.py",
    "scripts/render_ttf_q_figures.py",
    "tests/test_relational_qualification.py",
    "tests/test_relational_benchmark.py",
    "tests/test_ttf_q_known_truth_detectability_v02.py",
    "tests/test_ttf_q_bc_calibrated_detectability.py",
    "tests/test_ttf_q_manuscript_scalar_handoff.py",
    "tests/test_ttf_q_palearctic_prospective_stop.py",
]

PROHIBITED_TEXT = (
    "github.com/zuizui0223",
    "zuizui0223/TTF",
)

PROHIBITED_PATTERNS = {
    "public_github_url": re.compile(r"https?://(?:www\.)?github\.com/", re.I),
    "email_address": re.compile(
        r"\b[A-Z0-9._%+-]+@[A-Z0-9.-]+\.[A-Z]{2,}\b",
        re.I,
    ),
    "orcid_url": re.compile(r"https?://(?:www\.)?orcid\.org/", re.I),
}


def anonymity_hits(body: str) -> list[str]:
    hits = [
        f"literal:{token}"
        for token in PROHIBITED_TEXT
        if token.lower() in body.lower()
    ]
    hits.extend(
        f"pattern:{name}"
        for name, pattern in PROHIBITED_PATTERNS.items()
        if pattern.search(body)
    )
    return hits


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()

    license_path = ROOT / "LICENSE"
    if not license_path.is_file():
        raise RuntimeError(
            "LICENSE is required before building the MEE peer-review code archive. "
            "Select an open-source license explicitly; do not infer one."
        )

    missing = [name for name in FILES if not (ROOT / name).is_file()]
    if missing:
        raise RuntimeError(f"peer-review bundle inputs missing: {missing}")

    with tempfile.TemporaryDirectory() as td:
        stage = Path(td) / "ttf_q_anonymous_peer_review"
        stage.mkdir()
        for name in FILES:
            src = ROOT / name
            dst = stage / name
            dst.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(src, dst)

        pyproject = stage / "pyproject.toml"
        pyproject_text = pyproject.read_text()
        pyproject_text = pyproject_text.replace(
            'name = "transferable-turnover-fields"',
            'name = "ttf-q-peer-review"',
        )
        pyproject_text = pyproject_text.replace(
            'description = "Species-disjoint inference for transferable within-species trait-transition fields"',
            'description = "Anonymous peer-review implementation of TTF-Q"',
        )
        pyproject.write_text(pyproject_text)

        readme = stage / "README.md"
        readme.write_text(
            "# TTF-Q anonymous peer-review archive\n\n"
            "This archive contains the response-blind TTF-Q implementation, "
            "known-truth benchmark contracts and receipts, tests, manuscript "
            "scalar handoff, and figure renderer. Version-control history and "
            "public-repository metadata are intentionally excluded for "
            "double-anonymous review.\n\n"
            "Install with: python -m pip install -e .[test]\n\n"
            "Run the included focused tests with pytest.\n"
        )

        for path in stage.rglob("*"):
            if not path.is_file() or path.suffix.lower() in {".png", ".jpg", ".jpeg", ".zip"}:
                continue
            body = path.read_text(errors="ignore")
            hits = anonymity_hits(body)
            if hits:
                raise RuntimeError(
                    f"author-identifying public-repository token in "
                    f"{path.relative_to(stage)}: {hits}"
                )

        manifest = {
            "schema": "ttf_q_anonymous_peer_review_bundle_v0.1",
            "files": {
                str(path.relative_to(stage)): sha256(path)
                for path in sorted(stage.rglob("*"))
                if path.is_file()
            },
            "git_history_included": False,
            "public_repository_url_included": False,
            "license_included": True,
        }
        (stage / "BUNDLE_MANIFEST.json").write_text(
            json.dumps(manifest, indent=2, sort_keys=True) + "\n"
        )

        args.output.parent.mkdir(parents=True, exist_ok=True)
        with zipfile.ZipFile(args.output, "w", compression=zipfile.ZIP_DEFLATED) as zf:
            for path in sorted(stage.rglob("*")):
                if path.is_file():
                    zf.write(path, path.relative_to(stage.parent))

    print(json.dumps({
        "status": "BUILT_ANONYMOUS_PEER_REVIEW_BUNDLE",
        "output": str(args.output),
    }, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
