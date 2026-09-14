#!/usr/bin/env python3
from __future__ import annotations

import argparse
from contextlib import contextmanager
import json
import stat
import subprocess
import sys
import tarfile
import tempfile
from pathlib import Path
from typing import Iterator
import zipfile

from ttf.phylogatr_confirmatory import sha256_path


PROTOCOL = Path("docs/supporting/genetic_phylogatr_confirmatory_protocol_v0.1.json")
PARSER_RULE = Path("docs/supporting/genetic_phylogatr_phase1_parser_rule_v0.1.json")
DIGEST_RULE = Path("docs/supporting/genetic_phylogatr_phase1_digest_rule_v0.1.json")
EXECUTION_RULE = Path("docs/supporting/genetic_phylogatr_phase1_execution_rule_v0.1.json")
DECKER_EXCLUSION = Path("benchmarks/frozen/genetic_decker_species_exclusion_v0.1.json")
PHASE1_SCRIPT = Path("scripts/freeze_phylogatr_confirmatory_phase1.py")


def _repo_file(repo_root: Path, relative: Path) -> Path:
    path = repo_root / relative
    if not path.is_file():
        raise FileNotFoundError(f"required frozen Phase-1 file missing: {relative}")
    return path


def _safe_member_path(name: str) -> Path:
    path = Path(name)
    if path.is_absolute() or ".." in path.parts:
        raise RuntimeError(f"unsafe archive member path: {name!r}")
    return path


def _extract_zip(archive: Path, destination: Path) -> None:
    with zipfile.ZipFile(archive) as handle:
        for info in handle.infolist():
            _safe_member_path(info.filename)
            mode = (int(info.external_attr) >> 16) & 0o170000
            if mode == stat.S_IFLNK:
                raise RuntimeError(f"symlink archive member is forbidden: {info.filename!r}")
        handle.extractall(destination)


def _extract_tar(archive: Path, destination: Path) -> None:
    with tarfile.open(archive, mode="r:*") as handle:
        members = handle.getmembers()
        for member in members:
            _safe_member_path(member.name)
            if not (member.isfile() or member.isdir()):
                raise RuntimeError(
                    f"non-file/non-directory tar member is forbidden: {member.name!r}"
                )
        handle.extractall(destination, members=members, filter="data")


def _archive_format(path: Path) -> str:
    lower = path.name.lower()
    if lower.endswith(".zip"):
        return "zip"
    if lower.endswith(".tar.gz") or lower.endswith(".tgz"):
        return "tar.gz"
    raise ValueError("--archive must end in .zip, .tar.gz, or .tgz")


def _find_phylogatr_root(extracted: Path) -> Path:
    candidates: set[Path] = set()
    if (extracted / "genes.txt").is_file() and (extracted / "cite.txt").is_file():
        candidates.add(extracted)
    for genes in extracted.rglob("genes.txt"):
        parent = genes.parent
        if (parent / "cite.txt").is_file():
            candidates.add(parent)
    if len(candidates) != 1:
        shown = sorted(str(path.relative_to(extracted)) for path in candidates)
        raise RuntimeError(
            "archive must contain exactly one phylogatr-results root with genes.txt and cite.txt; "
            f"found {len(candidates)}: {shown[:10]}"
        )
    return next(iter(candidates))


@contextmanager
def _resolved_input_root(
    *,
    root: Path | None,
    archive: Path | None,
) -> Iterator[tuple[Path, dict[str, object]]]:
    if root is not None:
        resolved = root.resolve()
        if not (resolved / "genes.txt").is_file() or not (resolved / "cite.txt").is_file():
            raise FileNotFoundError(
                "--root must be a fresh phylogatr-results directory containing genes.txt and cite.txt"
            )
        yield resolved, {
            "source_type": "directory",
            "archive_name": None,
            "archive_sha256": None,
            "archive_format": None,
            "temporary_extraction": False,
        }
        return

    assert archive is not None
    resolved_archive = archive.resolve()
    if not resolved_archive.is_file():
        raise FileNotFoundError(f"fresh phylogatR archive not found: {resolved_archive}")
    archive_format = _archive_format(resolved_archive)
    archive_sha256 = sha256_path(resolved_archive)
    with tempfile.TemporaryDirectory(prefix="ttf_phylogatr_phase1_") as temporary:
        extracted = Path(temporary)
        if archive_format == "zip":
            _extract_zip(resolved_archive, extracted)
        else:
            _extract_tar(resolved_archive, extracted)
        discovered = _find_phylogatr_root(extracted)
        root_relative = str(discovered.relative_to(extracted)) or "."
        yield discovered, {
            "source_type": "archive",
            "archive_name": resolved_archive.name,
            "archive_sha256": archive_sha256,
            "archive_format": archive_format,
            "archive_root_relative": root_relative,
            "temporary_extraction": True,
            "temporary_extraction_deleted_after_phase1": True,
        }


def _run_phase1(
    *,
    root: Path,
    output_dir: Path,
    repo_root: Path,
    source_provenance: dict[str, object],
) -> int:
    protocol = _repo_file(repo_root, PROTOCOL)
    parser_rule = _repo_file(repo_root, PARSER_RULE)
    digest_rule = _repo_file(repo_root, DIGEST_RULE)
    execution_rule = _repo_file(repo_root, EXECUTION_RULE)
    exclusion = _repo_file(repo_root, DECKER_EXCLUSION)
    phase1_script = _repo_file(repo_root, PHASE1_SCRIPT)

    geometry_csv = output_dir / "phase1_geometry.csv"
    manifest_path = output_dir / "phase1_manifest.json"
    receipt_path = output_dir / "phase1_intake_receipt.json"

    if any(path.exists() for path in (geometry_csv, manifest_path, receipt_path)):
        raise FileExistsError(
            "Phase-1 intake outputs already exist; use a new output directory rather than overwriting a frozen receipt"
        )

    command = [
        sys.executable,
        str(phase1_script),
        "--root",
        str(root),
        "--protocol",
        str(protocol),
        "--parser-rule",
        str(parser_rule),
        "--digest-rule",
        str(digest_rule),
        "--execution-rule",
        str(execution_rule),
        "--decker-exclusion",
        str(exclusion),
        "--output-csv",
        str(geometry_csv),
        "--output-manifest",
        str(manifest_path),
    ]
    completed = subprocess.run(
        command,
        cwd=repo_root,
        check=False,
        text=True,
        capture_output=True,
    )
    if completed.returncode != 0:
        if geometry_csv.exists() or manifest_path.exists():
            raise RuntimeError(
                "Phase-1 freezer failed after creating partial outputs; do not reuse this output directory.\n"
                + completed.stderr
            )
        raise RuntimeError("Phase-1 freezer failed before freeze completion.\n" + completed.stderr)

    manifest = json.loads(manifest_path.read_text())
    if manifest.get("schema") != "ttf_genetic_phylogatr_confirmatory_phase1_geometry_v0.1":
        raise RuntimeError("unexpected Phase-1 manifest schema")
    if manifest.get("confirmatory_sequence_identity_opened") is not False:
        raise RuntimeError("Phase-1 manifest indicates nucleotide identity opening")
    if manifest.get("confirmatory_pairwise_genetic_distances_opened") is not False:
        raise RuntimeError("Phase-1 manifest indicates genetic-distance opening")
    response_blind = manifest.get("response_blind")
    if not isinstance(response_blind, dict):
        raise RuntimeError("Phase-1 response-blind receipt missing")
    if response_blind.get("sequence_characters_used") is not False:
        raise RuntimeError("Phase-1 intake used sequence characters")
    if response_blind.get("sequence_characters_hashed") is not False:
        raise RuntimeError("Phase-1 intake hashed sequence characters")

    passed = manifest.get("status") == "FROZEN_RESPONSE_BLIND_PHASE1_GEOMETRY"
    receipt = {
        "schema": "ttf_genetic_phylogatr_phase1_intake_receipt_v0.2",
        "status": "PHASE1_FROZEN" if passed else "NOT_EVALUABLE_PHASE1",
        "purpose": "Single-command receipt for the only currently authorized fresh-phylogatR action: response-blind Phase 1. Directory and safely extracted archive inputs are equivalent at the Phase-1 estimand surface.",
        "source_provenance": {
            **source_provenance,
            "cite_sha256": sha256_path(root / "cite.txt"),
            "genes_sha256": sha256_path(root / "genes.txt"),
        },
        "frozen_inputs_sha256": {
            str(PROTOCOL): sha256_path(protocol),
            str(PARSER_RULE): sha256_path(parser_rule),
            str(DIGEST_RULE): sha256_path(digest_rule),
            str(EXECUTION_RULE): sha256_path(execution_rule),
            str(DECKER_EXCLUSION): sha256_path(exclusion),
            str(PHASE1_SCRIPT): sha256_path(phase1_script),
        },
        "phase1_manifest_sha256": sha256_path(manifest_path),
        "phase1_geometry_csv_sha256": sha256_path(geometry_csv),
        "dataset_digest_sha256": manifest.get("dataset_digest_sha256"),
        "geometry_fingerprint_sha256": manifest.get("geometry_fingerprint_sha256"),
        "species": int(manifest.get("census", {}).get("final_species", 0)),
        "train_species": int(manifest.get("split", {}).get("train_count", 0)),
        "eval_species": int(manifest.get("split", {}).get("eval_count", 0)),
        "localities_total": int(manifest.get("localities_total", 0)),
        "edges_total": int(manifest.get("edges_total", 0)),
        "sequence_characters_opened": False,
        "sequence_characters_hashed": False,
        "pairwise_genetic_distances_opened": False,
        "empirical_ttf_opened": False,
        "next_authorized_action": (
            "Run the frozen Phase-2 canonical-valid/noncanonical mask admissibility on these exact Phase-1 files."
            if passed
            else "STOP. Do not open Phase-2 character masks or nucleotide identity."
        ),
        "phase2_automatic_execution": False,
    }
    receipt_path.write_text(json.dumps(receipt, indent=2, sort_keys=True) + "\n")

    print(
        json.dumps(
            {
                "status": receipt["status"],
                "source_type": receipt["source_provenance"]["source_type"],
                "species": receipt["species"],
                "train_species": receipt["train_species"],
                "eval_species": receipt["eval_species"],
                "dataset_digest_sha256": receipt["dataset_digest_sha256"],
                "geometry_fingerprint_sha256": receipt["geometry_fingerprint_sha256"],
                "sequence_characters_opened": False,
                "phase2_automatic_execution": False,
            },
            sort_keys=True,
        )
    )
    return 0


def main() -> int:
    ap = argparse.ArgumentParser(
        description=(
            "Run only the authorized response-blind Phase-1 intake on a fresh phylogatR directory or archive. "
            "Archive extraction copies bytes only; this wrapper never interprets nucleotide identity or opens Phase-2 masks."
        )
    )
    source = ap.add_mutually_exclusive_group(required=True)
    source.add_argument("--root", type=Path, help="fresh phylogatr-results directory")
    source.add_argument("--archive", type=Path, help="fresh phylogatR .zip, .tar.gz, or .tgz download")
    ap.add_argument("--output-dir", type=Path, required=True)
    ap.add_argument("--repo-root", type=Path, default=Path("."))
    args = ap.parse_args()

    repo_root = args.repo_root.resolve()
    output_dir = args.output_dir.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    with _resolved_input_root(root=args.root, archive=args.archive) as (root, provenance):
        return _run_phase1(
            root=root,
            output_dir=output_dir,
            repo_root=repo_root,
            source_provenance=provenance,
        )


if __name__ == "__main__":
    raise SystemExit(main())
