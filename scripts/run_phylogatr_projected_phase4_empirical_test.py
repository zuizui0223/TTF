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
INNER_PHASE4 = Path("scripts/run_phylogatr_phase4_empirical_test.py")


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
                    raise RuntimeError(
                        f"symlink archive member is forbidden: {info.filename!r}"
                    )
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
            "archive must contain exactly one phylogatR root with genes.txt and cite.txt; "
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


def _load_phase1_projection(phase1_dir: Path) -> tuple[Path, dict, dict]:
    manifest_path = phase1_dir / "phase1_manifest.json"
    receipt_path = phase1_dir / "phase1_intake_receipt.json"
    if not manifest_path.is_file() or not receipt_path.is_file():
        raise FileNotFoundError("projected Phase-1 manifest/intake receipt is missing")

    manifest = json.loads(manifest_path.read_text())
    if manifest.get("schema") != "ttf_genetic_phylogatr_confirmatory_phase1_geometry_v0.1":
        raise RuntimeError("unexpected Phase-1 manifest schema")

    receipt = json.loads(receipt_path.read_text())
    if receipt.get("schema") != "ttf_genetic_phylogatr_phase1_intake_receipt_v0.2":
        raise RuntimeError("unexpected Phase-1 intake receipt schema")
    if receipt.get("status") != "PHASE1_FROZEN":
        raise RuntimeError("projected Phase-4 requires a passed Phase-1 receipt")

    phase1_manifest_sha256 = sha256_path(manifest_path)
    if receipt.get("phase1_manifest_sha256") != phase1_manifest_sha256:
        raise RuntimeError("Phase-1 intake receipt/manifest SHA256 linkage drift")

    projection = receipt.get("source_projection")
    if not isinstance(projection, dict):
        raise RuntimeError("Phase-1 receipt lacks source_projection provenance")
    if projection.get("schema") != "ttf_genetic_phylogatr_source_projection_receipt_v0.1":
        raise RuntimeError("unexpected Phase-1 source_projection schema")
    if projection.get("status") != "ANIMALIA_METADATA_PROJECTION_APPLIED":
        raise RuntimeError("Phase-1 source_projection was not frozen as Animalia projection")
    if projection.get("sequence_identity_opened") is not False:
        raise RuntimeError("Phase-1 source_projection indicates identity opening")
    if manifest.get("provenance", {}).get("genes_sha256") != projection.get(
        "projected_genes_sha256"
    ):
        raise RuntimeError("Phase-1 manifest/projection genes SHA256 drift")
    if manifest.get("provenance", {}).get("cite_sha256") != projection.get("cite_sha256"):
        raise RuntimeError("Phase-1 manifest/projection cite SHA256 drift")
    return manifest_path, receipt, projection


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
            raise RuntimeError(f"Phase-4 source projection differs from Phase 1: {key}")


def main() -> int:
    ap = argparse.ArgumentParser(
        description=(
            "Run the one authorized fresh phylogatR empirical test from the exact broader "
            "authenticated archive by recreating the frozen pre-identity Animalia metadata "
            "projection in a temporary source view before nucleotide identity is opened."
        )
    )
    ap.add_argument("--archive", type=Path, required=True)
    ap.add_argument("--phase1-dir", type=Path, required=True)
    ap.add_argument("--geometry", type=Path, required=True)
    ap.add_argument("--phase2-manifest", type=Path, required=True)
    ap.add_argument("--phase3-rule", type=Path, required=True)
    ap.add_argument("--phase3-authorization", type=Path, required=True)
    ap.add_argument("--references", type=Path, required=True)
    ap.add_argument("--qualification", type=Path, required=True)
    ap.add_argument("--self-rule", type=Path, required=True)
    ap.add_argument("--self-references", type=Path, required=True)
    ap.add_argument("--self-qualification", type=Path, required=True)
    ap.add_argument("--phase4-rule", type=Path, required=True)
    ap.add_argument("--phase4-authorization", type=Path, required=True)
    ap.add_argument("--opening-state", type=Path, required=True)
    ap.add_argument("--repo-root", type=Path, default=Path("."))
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()

    repo_root = args.repo_root.resolve()
    archive = args.archive.resolve()
    phase1_dir = args.phase1_dir.resolve()
    output = args.output.resolve()
    if not archive.is_file():
        raise FileNotFoundError(f"phylogatR archive not found: {archive}")

    projection_rule = repo_root / PROJECTION_RULE
    inner_phase4 = repo_root / INNER_PHASE4
    if not projection_rule.is_file() or not inner_phase4.is_file():
        raise FileNotFoundError("required projected Phase-4 frozen files are missing")
    _load_projection_rule(projection_rule)
    phase1_manifest, phase1_receipt, frozen_projection = _load_phase1_projection(
        phase1_dir
    )

    raw_archive_sha256 = sha256_path(archive)
    if raw_archive_sha256 != frozen_projection.get("raw_archive_sha256"):
        raise RuntimeError("Phase-4 raw archive SHA256 differs from projected Phase 1")
    archive_format = _archive_format(archive)
    if archive_format != frozen_projection.get("raw_archive_format"):
        raise RuntimeError("Phase-4 raw archive format differs from projected Phase 1")
    projection_rule_sha256 = sha256_path(projection_rule)
    if projection_rule_sha256 != frozen_projection.get("projection_rule_sha256"):
        raise RuntimeError("source-projection rule SHA256 drift since Phase 1")

    empirical_result: dict
    projection: dict
    with tempfile.TemporaryDirectory(prefix="ttf_phylogatr_projected_phase4_") as temporary:
        extracted = Path(temporary)
        _extract_archive(archive, extracted, archive_format)
        root = _find_phylogatr_root(extracted)
        root_relative = str(root.relative_to(extracted)) or "."
        if root_relative != frozen_projection.get("archive_root_relative"):
            raise RuntimeError("Phase-4 archive root differs from projected Phase 1")

        projection = project_genes_file_in_place(root)
        _require_same_projection(projection, frozen_projection)
        if projection["projected_genes_sha256"] != json.loads(
            phase1_manifest.read_text()
        )["provenance"]["genes_sha256"]:
            raise RuntimeError("Phase-4 projected genes.txt does not match Phase-1 manifest")
        projection.update(
            {
                "raw_archive_name": archive.name,
                "raw_archive_sha256": raw_archive_sha256,
                "raw_archive_format": archive_format,
                "archive_root_relative": root_relative,
                "projection_rule": str(PROJECTION_RULE),
                "projection_rule_sha256": projection_rule_sha256,
                "phase1_manifest_sha256": sha256_path(phase1_manifest),
                "phase1_intake_receipt_sha256": sha256_path(
                    phase1_dir / "phase1_intake_receipt.json"
                ),
                "temporary_projection_view": True,
                "temporary_projection_view_deleted_after_phase4": True,
            }
        )

        command = [
            sys.executable,
            str(inner_phase4),
            "--root",
            str(root),
            "--geometry",
            str(args.geometry.resolve()),
            "--phase1-manifest",
            str(phase1_manifest),
            "--phase2-manifest",
            str(args.phase2_manifest.resolve()),
            "--phase3-rule",
            str(args.phase3_rule.resolve()),
            "--phase3-authorization",
            str(args.phase3_authorization.resolve()),
            "--references",
            str(args.references.resolve()),
            "--qualification",
            str(args.qualification.resolve()),
            "--self-rule",
            str(args.self_rule.resolve()),
            "--self-references",
            str(args.self_references.resolve()),
            "--self-qualification",
            str(args.self_qualification.resolve()),
            "--phase4-rule",
            str(args.phase4_rule.resolve()),
            "--phase4-authorization",
            str(args.phase4_authorization.resolve()),
            "--opening-state",
            str(args.opening_state.resolve()),
            "--repo-root",
            str(repo_root),
            "--output",
            str(output),
        ]
        completed = subprocess.run(
            command,
            cwd=repo_root,
            check=False,
            text=True,
            capture_output=True,
        )
        if completed.returncode != 0:
            raise RuntimeError(
                "projected Phase-4 inner empirical test failed.\n"
                + completed.stdout
                + "\n"
                + completed.stderr
            )
        if not output.is_file():
            raise RuntimeError("projected Phase-4 inner empirical test wrote no result")
        empirical_result = json.loads(output.read_text())
        if (
            empirical_result.get("schema")
            != "ttf_genetic_phylogatr_phase4_empirical_result_v0.1"
        ):
            raise RuntimeError("unexpected Phase-4 empirical result schema")

    if empirical_result.get("confirmatory_sequence_identity_opened") is not True:
        raise RuntimeError("Phase-4 empirical result did not record identity opening")
    if empirical_result.get("confirmatory_pairwise_genetic_distances_opened") is not True:
        raise RuntimeError("Phase-4 empirical result did not record distance opening")
    if empirical_result.get("confirmatory_ttf_statistic_opened") is not True:
        raise RuntimeError("Phase-4 empirical result did not record TTF opening")

    empirical_result["source_projection"] = projection
    empirical_result["source_projection"]["sequence_identity_opened"] = True
    empirical_result["source_projection"]["pairwise_genetic_distances_opened"] = True
    empirical_result["source_projection"]["empirical_ttf_opened"] = True
    empirical_result["source_projection"]["projection_itself_opened_identity"] = False
    output.write_text(
        json.dumps(empirical_result, indent=2, sort_keys=True, allow_nan=False) + "\n"
    )

    print(
        json.dumps(
            {
                "status": empirical_result["status"],
                "decision": empirical_result["decision"],
                "raw_archive_sha256": raw_archive_sha256,
                "projected_genes_sha256": projection["projected_genes_sha256"],
                "phase1_manifest_sha256": projection["phase1_manifest_sha256"],
                "sequence_identity_opened": True,
                "pairwise_genetic_distances_opened": True,
                "empirical_ttf_opened": True,
            },
            sort_keys=True,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
