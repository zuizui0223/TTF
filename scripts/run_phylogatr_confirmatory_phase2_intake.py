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


PHASE2_RULE = Path("docs/supporting/genetic_phylogatr_phase2_mask_rule_v0.1.json")
PHASE2_SCRIPT = Path("scripts/freeze_phylogatr_confirmatory_phase2_mask.py")


def _repo_file(repo_root: Path, relative: Path) -> Path:
    path = repo_root / relative
    if not path.is_file():
        raise FileNotFoundError(f"required frozen Phase-2 file missing: {relative}")
    return path


def _load_json(path: Path, schema: str) -> dict:
    payload = json.loads(path.read_text())
    if payload.get("schema") != schema:
        raise RuntimeError(f"unexpected schema for {path}: {payload.get('schema')!r}")
    return payload


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


def _validate_phase1_chain(phase1_dir: Path) -> tuple[dict, dict, Path, Path, Path]:
    manifest_path = phase1_dir / "phase1_manifest.json"
    geometry_path = phase1_dir / "phase1_geometry.csv"
    receipt_path = phase1_dir / "phase1_intake_receipt.json"
    for path in (manifest_path, geometry_path, receipt_path):
        if not path.is_file():
            raise FileNotFoundError(f"required frozen Phase-1 output missing: {path}")

    manifest = _load_json(
        manifest_path,
        "ttf_genetic_phylogatr_confirmatory_phase1_geometry_v0.1",
    )
    receipt = _load_json(
        receipt_path,
        "ttf_genetic_phylogatr_phase1_intake_receipt_v0.2",
    )
    if manifest.get("status") != "FROZEN_RESPONSE_BLIND_PHASE1_GEOMETRY":
        raise RuntimeError("Phase-2 intake requires a passed frozen Phase-1 geometry")
    if receipt.get("status") != "PHASE1_FROZEN":
        raise RuntimeError("Phase-2 intake requires a PHASE1_FROZEN intake receipt")
    if sha256_path(manifest_path) != receipt.get("phase1_manifest_sha256"):
        raise RuntimeError("Phase-1 manifest SHA256 drift since intake")
    if sha256_path(geometry_path) != receipt.get("phase1_geometry_csv_sha256"):
        raise RuntimeError("Phase-1 geometry CSV SHA256 drift since intake")
    if receipt.get("dataset_digest_sha256") != manifest.get("dataset_digest_sha256"):
        raise RuntimeError("Phase-1 dataset digest drift between manifest and intake receipt")
    if receipt.get("geometry_fingerprint_sha256") != manifest.get("geometry_fingerprint_sha256"):
        raise RuntimeError("Phase-1 geometry fingerprint drift between manifest and intake receipt")
    if manifest.get("confirmatory_sequence_identity_opened") is not False:
        raise RuntimeError("Phase-1 manifest indicates nucleotide identity was opened")
    if manifest.get("confirmatory_pairwise_genetic_distances_opened") is not False:
        raise RuntimeError("Phase-1 manifest indicates genetic distances were opened")
    return manifest, receipt, manifest_path, geometry_path, receipt_path


@contextmanager
def _resolved_input_root(
    *,
    root: Path | None,
    archive: Path | None,
    phase1_receipt: dict,
) -> Iterator[tuple[Path, dict[str, object]]]:
    expected = phase1_receipt.get("source_provenance")
    if not isinstance(expected, dict):
        raise RuntimeError("Phase-1 intake receipt lacks source provenance")
    expected_type = expected.get("source_type")

    if root is not None:
        if expected_type != "directory":
            raise RuntimeError(
                "Phase-1 used an archive; Phase-2 must receive the same untouched archive, not a manually extracted directory"
            )
        resolved = root.resolve()
        if not (resolved / "genes.txt").is_file() or not (resolved / "cite.txt").is_file():
            raise FileNotFoundError(
                "--root must be the same fresh phylogatr-results directory used for Phase 1"
            )
        yield resolved, {
            "source_type": "directory",
            "temporary_extraction": False,
            "phase1_source_type_matched": True,
        }
        return

    if expected_type != "archive":
        raise RuntimeError(
            "Phase-1 used a directory; Phase-2 must receive that same directory rather than a new archive"
        )
    assert archive is not None
    resolved_archive = archive.resolve()
    if not resolved_archive.is_file():
        raise FileNotFoundError(f"fresh phylogatR archive not found: {resolved_archive}")
    archive_format = _archive_format(resolved_archive)
    archive_sha256 = sha256_path(resolved_archive)
    if archive_sha256 != expected.get("archive_sha256"):
        raise RuntimeError("Phase-2 archive SHA256 differs from the untouched Phase-1 archive")
    if archive_format != expected.get("archive_format"):
        raise RuntimeError("Phase-2 archive format differs from Phase-1 provenance")

    with tempfile.TemporaryDirectory(prefix="ttf_phylogatr_phase2_") as temporary:
        extracted = Path(temporary)
        if archive_format == "zip":
            _extract_zip(resolved_archive, extracted)
        else:
            _extract_tar(resolved_archive, extracted)
        discovered = _find_phylogatr_root(extracted)
        root_relative = str(discovered.relative_to(extracted)) or "."
        if root_relative != expected.get("archive_root_relative"):
            raise RuntimeError("Phase-2 archive root differs from Phase-1 archive provenance")
        yield discovered, {
            "source_type": "archive",
            "archive_name": resolved_archive.name,
            "archive_sha256": archive_sha256,
            "archive_format": archive_format,
            "archive_root_relative": root_relative,
            "temporary_extraction": True,
            "temporary_extraction_deleted_after_phase2": True,
            "phase1_source_type_matched": True,
        }


def _run_phase2(
    *,
    root: Path,
    phase1_manifest_path: Path,
    phase1_geometry_path: Path,
    phase1_receipt_path: Path,
    phase1_manifest: dict,
    phase1_receipt: dict,
    output_dir: Path,
    repo_root: Path,
    source_provenance: dict[str, object],
) -> int:
    phase2_rule = _repo_file(repo_root, PHASE2_RULE)
    phase2_script = _repo_file(repo_root, PHASE2_SCRIPT)
    output_csv = output_dir / "phase2_geometry.csv"
    output_manifest = output_dir / "phase2_manifest.json"
    output_receipt = output_dir / "phase2_intake_receipt.json"
    if any(path.exists() for path in (output_csv, output_manifest, output_receipt)):
        raise FileExistsError(
            "Phase-2 intake outputs already exist; use a new output directory rather than overwriting a frozen receipt"
        )

    expected_source = phase1_receipt["source_provenance"]
    if sha256_path(root / "cite.txt") != expected_source.get("cite_sha256"):
        raise RuntimeError("cite.txt differs from Phase-1 source provenance")
    if sha256_path(root / "genes.txt") != expected_source.get("genes_sha256"):
        raise RuntimeError("genes.txt differs from Phase-1 source provenance")

    command = [
        sys.executable,
        str(phase2_script),
        "--root",
        str(root),
        "--phase1-manifest",
        str(phase1_manifest_path),
        "--phase1-csv",
        str(phase1_geometry_path),
        "--phase2-rule",
        str(phase2_rule),
        "--output-csv",
        str(output_csv),
        "--output-manifest",
        str(output_manifest),
    ]
    completed = subprocess.run(
        command,
        cwd=repo_root,
        check=False,
        text=True,
        capture_output=True,
    )
    if completed.returncode != 0:
        if output_csv.exists() or output_manifest.exists():
            raise RuntimeError(
                "Phase-2 freezer failed after creating partial outputs; do not reuse this output directory.\n"
                + completed.stderr
            )
        raise RuntimeError("Phase-2 freezer failed before freeze completion.\n" + completed.stderr)

    manifest = _load_json(
        output_manifest,
        "ttf_genetic_phylogatr_confirmatory_phase2_mask_v0.1",
    )
    if manifest.get("character_mask_opened") is not True:
        raise RuntimeError("Phase-2 manifest does not record character-mask opening")
    if manifest.get("confirmatory_sequence_identity_opened") is not False:
        raise RuntimeError("Phase-2 manifest indicates nucleotide identity opening")
    if manifest.get("confirmatory_pairwise_genetic_distances_opened") is not False:
        raise RuntimeError("Phase-2 manifest indicates pairwise genetic-distance opening")
    if manifest.get("confirmatory_ttf_statistic_opened") is not False:
        raise RuntimeError("Phase-2 manifest indicates empirical TTF opening")
    passed = manifest.get("status") == "PASS_TO_SYNTHETIC_GATE"
    receipt = {
        "schema": "ttf_genetic_phylogatr_phase2_intake_receipt_v0.1",
        "status": "PHASE2_FROZEN_TO_SYNTHETIC_GATE" if passed else "NOT_EVALUABLE_PHASE2",
        "purpose": "Archive-safe Phase-2 canonical-valid/noncanonical mask intake bound to the exact Phase-1 source. The wrapper may reveal only the frozen character-validity mask; nucleotide identity remains closed.",
        "source_provenance": {
            **source_provenance,
            "cite_sha256": sha256_path(root / "cite.txt"),
            "genes_sha256": sha256_path(root / "genes.txt"),
        },
        "phase1_intake_receipt_sha256": sha256_path(phase1_receipt_path),
        "phase1_manifest_sha256": sha256_path(phase1_manifest_path),
        "phase1_geometry_csv_sha256": sha256_path(phase1_geometry_path),
        "phase2_rule_sha256": sha256_path(phase2_rule),
        "phase2_script_sha256": sha256_path(phase2_script),
        "phase2_manifest_sha256": sha256_path(output_manifest),
        "phase2_geometry_csv_sha256": sha256_path(output_csv),
        "dataset_digest_sha256": phase1_manifest["dataset_digest_sha256"],
        "phase1_geometry_fingerprint_sha256": phase1_manifest["geometry_fingerprint_sha256"],
        "phase2_geometry_fingerprint_sha256": manifest.get("geometry_fingerprint_sha256"),
        "phase1_species": int(manifest.get("species", {}).get("phase1", 0)),
        "survivor_species": int(manifest.get("species", {}).get("survivors", 0)),
        "failed_character_support": int(manifest.get("species", {}).get("failed_character_support", 0)),
        "train_species": int(manifest.get("split", {}).get("train_count", 0)),
        "eval_species": int(manifest.get("split", {}).get("eval_count", 0)),
        "character_mask_opened": True,
        "sequence_identity_opened": False,
        "pairwise_genetic_distances_opened": False,
        "empirical_ttf_opened": False,
        "next_authorized_action": (
            "Authorize and run the frozen formal Phase-3 Gate-D on this exact Phase-2 survivor geometry."
            if passed
            else "STOP. Do not open nucleotide identity or pairwise genetic distances."
        ),
        "phase3_automatic_execution": False,
    }
    output_receipt.write_text(json.dumps(receipt, indent=2, sort_keys=True) + "\n")
    print(
        json.dumps(
            {
                "status": receipt["status"],
                "source_type": receipt["source_provenance"]["source_type"],
                "survivor_species": receipt["survivor_species"],
                "train_species": receipt["train_species"],
                "eval_species": receipt["eval_species"],
                "character_mask_opened": True,
                "sequence_identity_opened": False,
                "phase3_automatic_execution": False,
            },
            sort_keys=True,
        )
    )
    return 0


def main() -> int:
    ap = argparse.ArgumentParser(
        description=(
            "Run only the frozen Phase-2 character-mask intake on the exact source used by Phase 1. "
            "Archive inputs are re-extracted safely into a temporary directory and verified against the Phase-1 archive SHA256."
        )
    )
    source = ap.add_mutually_exclusive_group(required=True)
    source.add_argument("--root", type=Path, help="same fresh phylogatr-results directory used for Phase 1")
    source.add_argument("--archive", type=Path, help="same untouched .zip/.tar.gz/.tgz used for Phase 1")
    ap.add_argument("--phase1-dir", type=Path, required=True)
    ap.add_argument("--output-dir", type=Path, required=True)
    ap.add_argument("--repo-root", type=Path, default=Path("."))
    args = ap.parse_args()

    repo_root = args.repo_root.resolve()
    phase1_dir = args.phase1_dir.resolve()
    output_dir = args.output_dir.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    manifest, receipt, manifest_path, geometry_path, receipt_path = _validate_phase1_chain(
        phase1_dir
    )
    with _resolved_input_root(
        root=args.root,
        archive=args.archive,
        phase1_receipt=receipt,
    ) as (root, provenance):
        return _run_phase2(
            root=root,
            phase1_manifest_path=manifest_path,
            phase1_geometry_path=geometry_path,
            phase1_receipt_path=receipt_path,
            phase1_manifest=manifest,
            phase1_receipt=receipt,
            output_dir=output_dir,
            repo_root=repo_root,
            source_provenance=provenance,
        )


if __name__ == "__main__":
    raise SystemExit(main())
