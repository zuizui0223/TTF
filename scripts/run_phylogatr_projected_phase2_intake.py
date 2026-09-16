#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import stat
import subprocess
import sys
import tarfile
import tempfile
from pathlib import Path
import zipfile

from ttf.phylogatr_confirmatory import sha256_path
from ttf.phylogatr_source_projection import project_genes_file_in_place


PROJECTION_RULE = Path("docs/supporting/genetic_phylogatr_source_projection_rule_v0.1.json")
INNER_PHASE2 = Path("scripts/run_phylogatr_confirmatory_phase2_intake.py")


def _safe_member_path(name: str) -> Path:
    path = Path(name)
    if path.is_absolute() or ".." in path.parts:
        raise RuntimeError(f"unsafe archive member path: {name!r}")
    return path


def _archive_format(path: Path) -> str:
    lower = path.name.lower()
    if lower.endswith(".zip"):
        return "zip"
    if lower.endswith(".tar.gz") or lower.endswith(".tgz"):
        return "tar.gz"
    raise ValueError("--archive must end in .zip, .tar.gz, or .tgz")


def _extract_archive(archive: Path, destination: Path, archive_format: str) -> None:
    if archive_format == "zip":
        with zipfile.ZipFile(archive) as handle:
            for info in handle.infolist():
                _safe_member_path(info.filename)
                mode = (int(info.external_attr) >> 16) & 0o170000
                if mode == stat.S_IFLNK:
                    raise RuntimeError(f"symlink archive member is forbidden: {info.filename!r}")
            handle.extractall(destination)
        return
    with tarfile.open(archive, mode="r:*") as handle:
        members = handle.getmembers()
        for member in members:
            _safe_member_path(member.name)
            if not (member.isfile() or member.isdir()):
                raise RuntimeError(
                    f"non-file/non-directory tar member is forbidden: {member.name!r}"
                )
        handle.extractall(destination, members=members, filter="data")


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


def _load_projection_rule(path: Path) -> dict:
    payload = json.loads(path.read_text())
    if payload.get("schema") != "ttf_genetic_phylogatr_source_projection_rule_v0.1":
        raise RuntimeError("unexpected phylogatR source-projection rule schema")
    firewall = payload.get("firewall")
    if not isinstance(firewall, dict) or any(value is not False for value in firewall.values()):
        raise RuntimeError("source-projection outcome firewall is open")
    if payload.get("scientific_authority", {}).get("can_open_identity") is not False:
        raise RuntimeError("source-projection rule must not authorize identity opening")
    return payload


def _load_phase1_projection(phase1_dir: Path) -> tuple[dict, dict]:
    receipt_path = phase1_dir / "phase1_intake_receipt.json"
    if not receipt_path.is_file():
        raise FileNotFoundError(f"projected Phase-1 receipt missing: {receipt_path}")
    receipt = json.loads(receipt_path.read_text())
    if receipt.get("schema") != "ttf_genetic_phylogatr_phase1_intake_receipt_v0.2":
        raise RuntimeError("unexpected Phase-1 intake receipt schema")
    if receipt.get("status") != "PHASE1_FROZEN":
        raise RuntimeError("projected Phase-2 requires a passed Phase-1 receipt")
    projection = receipt.get("source_projection")
    if not isinstance(projection, dict):
        raise RuntimeError("Phase-1 receipt lacks source_projection provenance")
    if projection.get("schema") != "ttf_genetic_phylogatr_source_projection_receipt_v0.1":
        raise RuntimeError("unexpected Phase-1 source_projection schema")
    if projection.get("status") != "ANIMALIA_METADATA_PROJECTION_APPLIED":
        raise RuntimeError("Phase-1 source_projection was not frozen as Animalia projection")
    if projection.get("sequence_identity_opened") is not False:
        raise RuntimeError("Phase-1 source_projection indicates identity opening")
    return receipt, projection


def _require_same_projection(current: dict, frozen: dict) -> None:
    keys = (
        "raw_genes_sha256",
        "projected_genes_sha256",
        "cite_sha256",
        "original_row_count",
        "retained_animalia_row_count",
        "removed_non_animalia_row_count",
        "kingdom_counts",
    )
    for key in keys:
        if current.get(key) != frozen.get(key):
            raise RuntimeError(f"Phase-2 source projection differs from Phase 1: {key}")


def main() -> int:
    ap = argparse.ArgumentParser(
        description=(
            "Run the frozen Phase-2 character-mask intake on the exact authenticated archive and "
            "deterministic Animalia metadata projection used for projected Phase 1. Nucleotide "
            "identity and genetic distances remain closed."
        )
    )
    ap.add_argument("--archive", type=Path, required=True)
    ap.add_argument("--phase1-dir", type=Path, required=True)
    ap.add_argument("--output-dir", type=Path, required=True)
    ap.add_argument("--repo-root", type=Path, default=Path("."))
    args = ap.parse_args()

    repo_root = args.repo_root.resolve()
    archive = args.archive.resolve()
    phase1_dir = args.phase1_dir.resolve()
    output_dir = args.output_dir.resolve()
    if not archive.is_file():
        raise FileNotFoundError(f"phylogatR archive not found: {archive}")
    output_dir.mkdir(parents=True, exist_ok=True)

    projection_rule = repo_root / PROJECTION_RULE
    inner_phase2 = repo_root / INNER_PHASE2
    if not projection_rule.is_file() or not inner_phase2.is_file():
        raise FileNotFoundError("required projected Phase-2 frozen files are missing")
    _load_projection_rule(projection_rule)
    phase1_receipt, frozen_projection = _load_phase1_projection(phase1_dir)

    raw_archive_sha256 = sha256_path(archive)
    if raw_archive_sha256 != frozen_projection.get("raw_archive_sha256"):
        raise RuntimeError("Phase-2 raw archive SHA256 differs from projected Phase 1")
    archive_format = _archive_format(archive)
    if archive_format != frozen_projection.get("raw_archive_format"):
        raise RuntimeError("Phase-2 raw archive format differs from projected Phase 1")
    projection_rule_sha256 = sha256_path(projection_rule)
    if projection_rule_sha256 != frozen_projection.get("projection_rule_sha256"):
        raise RuntimeError("source-projection rule SHA256 drift since Phase 1")

    with tempfile.TemporaryDirectory(prefix="ttf_phylogatr_projected_phase2_") as temporary:
        extracted = Path(temporary)
        _extract_archive(archive, extracted, archive_format)
        root = _find_phylogatr_root(extracted)
        root_relative = str(root.relative_to(extracted)) or "."
        if root_relative != frozen_projection.get("archive_root_relative"):
            raise RuntimeError("Phase-2 archive root differs from projected Phase 1")

        projection = project_genes_file_in_place(root)
        _require_same_projection(projection, frozen_projection)
        projection.update(
            {
                "raw_archive_name": archive.name,
                "raw_archive_sha256": raw_archive_sha256,
                "raw_archive_format": archive_format,
                "archive_root_relative": root_relative,
                "projection_rule": str(PROJECTION_RULE),
                "projection_rule_sha256": projection_rule_sha256,
                "temporary_projection_view": True,
                "temporary_projection_view_deleted_after_phase2": True,
            }
        )

        command = [
            sys.executable,
            str(inner_phase2),
            "--root",
            str(root),
            "--phase1-dir",
            str(phase1_dir),
            "--output-dir",
            str(output_dir),
            "--repo-root",
            str(repo_root),
        ]
        completed = subprocess.run(
            command,
            cwd=repo_root,
            check=False,
            text=True,
            capture_output=True,
        )
        if completed.returncode != 0:
            raise RuntimeError("projected Phase-2 inner intake failed.\n" + completed.stderr)

        receipt_path = output_dir / "phase2_intake_receipt.json"
        if not receipt_path.is_file():
            raise RuntimeError("projected Phase-2 inner intake did not write its receipt")
        receipt = json.loads(receipt_path.read_text())
        if receipt.get("character_mask_opened") is not True:
            raise RuntimeError("projected Phase-2 did not record character-mask opening")
        if receipt.get("sequence_identity_opened") is not False:
            raise RuntimeError("projected Phase-2 opened nucleotide identity")
        if receipt.get("pairwise_genetic_distances_opened") is not False:
            raise RuntimeError("projected Phase-2 opened genetic distances")
        if receipt.get("empirical_ttf_opened") is not False:
            raise RuntimeError("projected Phase-2 opened empirical TTF")

        source = receipt.get("source_provenance")
        if not isinstance(source, dict):
            raise RuntimeError("Phase-2 receipt source provenance missing")
        if source.get("genes_sha256") != projection["projected_genes_sha256"]:
            raise RuntimeError("Phase-2 receipt did not use the projected genes.txt")
        if source.get("cite_sha256") != projection["cite_sha256"]:
            raise RuntimeError("Phase-2 cite.txt provenance drift after projection")

        receipt["source_projection"] = projection
        receipt["purpose"] = (
            str(receipt.get("purpose", ""))
            + " Source amendment: exact same authenticated archive and deterministic Animalia metadata projection as projected Phase 1."
        ).strip()
        receipt_path.write_text(json.dumps(receipt, indent=2, sort_keys=True) + "\n")

    print(
        json.dumps(
            {
                "status": receipt["status"],
                "raw_archive_sha256": raw_archive_sha256,
                "projected_genes_sha256": projection["projected_genes_sha256"],
                "survivor_species": receipt.get("survivor_species"),
                "character_mask_opened": True,
                "sequence_identity_opened": False,
            },
            sort_keys=True,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
