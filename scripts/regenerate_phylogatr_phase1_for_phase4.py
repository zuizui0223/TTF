#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import subprocess
import sys
from pathlib import Path

from ttf.phylogatr_confirmatory import sha256_path


REGENERATION_RULE = Path(
    "docs/supporting/genetic_phylogatr_phase1_source_regeneration_v0.1.json"
)
PROJECTED_PHASE1 = Path("scripts/run_phylogatr_projected_phase1_intake.py")


def _load_rule(path: Path) -> dict:
    payload = json.loads(path.read_text())
    if payload.get("schema") != "ttf_genetic_phylogatr_phase1_source_regeneration_v0.1":
        raise RuntimeError("unexpected Phase-1 source-regeneration rule schema")
    if payload.get("status") != "FROZEN_BEFORE_PHASE4_IDENTITY_OPENING":
        raise RuntimeError("Phase-1 source-regeneration rule is not frozen")
    firewall = payload.get("outcome_firewall")
    if not isinstance(firewall, dict) or not firewall or any(
        value is not False for value in firewall.values()
    ):
        raise RuntimeError("Phase-1 source-regeneration outcome firewall is open")
    authority = payload.get("authority")
    if not isinstance(authority, dict) or not authority or any(
        value is not False for value in authority.values()
    ):
        raise RuntimeError("Phase-1 source-regeneration rule gained authority")
    return payload


def validate_regenerated_phase1(
    archive: Path,
    output_dir: Path,
    rule: dict,
) -> dict:
    archive = Path(archive)
    output_dir = Path(output_dir)
    manifest_path = output_dir / "phase1_manifest.json"
    geometry_path = output_dir / "phase1_geometry.csv"
    intake_path = output_dir / "phase1_intake_receipt.json"
    for path in (manifest_path, geometry_path, intake_path):
        if not path.is_file():
            raise FileNotFoundError(f"regenerated Phase-1 output missing: {path.name}")

    raw = rule["raw_source"]
    if archive.stat().st_size != int(raw["archive_size_bytes"]):
        raise RuntimeError("raw phylogatR archive size differs from frozen source")
    raw_sha = sha256_path(archive)
    if raw_sha != str(raw["raw_archive_sha256"]):
        raise RuntimeError("raw phylogatR archive SHA256 differs from frozen source")

    expected = rule["expected_phase1"]
    manifest_sha = sha256_path(manifest_path)
    geometry_sha = sha256_path(geometry_path)
    if manifest_sha != str(expected["phase1_manifest_sha256"]):
        raise RuntimeError("regenerated Phase-1 manifest is not byte-identical to frozen Phase 1")
    if geometry_sha != str(expected["phase1_geometry_csv_sha256"]):
        raise RuntimeError("regenerated Phase-1 geometry is not byte-identical to frozen Phase 1")

    manifest = json.loads(manifest_path.read_text())
    if manifest.get("schema") != "ttf_genetic_phylogatr_confirmatory_phase1_geometry_v0.1":
        raise RuntimeError("unexpected regenerated Phase-1 manifest schema")
    if manifest.get("status") != "FROZEN_RESPONSE_BLIND_PHASE1_GEOMETRY":
        raise RuntimeError("regenerated Phase-1 panel did not pass its frozen geometry gate")
    if manifest.get("dataset_digest_sha256") != expected["dataset_digest_sha256"]:
        raise RuntimeError("regenerated Phase-1 dataset digest drift")
    if manifest.get("geometry_fingerprint_sha256") != expected["geometry_fingerprint_sha256"]:
        raise RuntimeError("regenerated Phase-1 geometry fingerprint drift")
    census = manifest.get("census", {})
    split = manifest.get("split", {})
    if int(census.get("final_species", -1)) != int(expected["species_count"]):
        raise RuntimeError("regenerated Phase-1 species count drift")
    if int(split.get("train_count", -1)) != int(expected["train_count"]):
        raise RuntimeError("regenerated Phase-1 training count drift")
    if int(split.get("eval_count", -1)) != int(expected["eval_count"]):
        raise RuntimeError("regenerated Phase-1 evaluation count drift")
    response_blind = manifest.get("response_blind", {})
    for key in (
        "sequence_characters_used",
        "sequence_characters_hashed",
        "pairwise_genetic_distances_opened",
        "ttf_statistic_opened",
        "decker_empirical_genetic_outcomes_opened",
    ):
        if response_blind.get(key) is not False:
            raise RuntimeError(f"regenerated Phase-1 response-blind firewall open: {key}")
    if manifest.get("confirmatory_sequence_identity_opened") is not False:
        raise RuntimeError("regenerated Phase-1 manifest indicates identity opening")
    if manifest.get("confirmatory_pairwise_genetic_distances_opened") is not False:
        raise RuntimeError("regenerated Phase-1 manifest indicates distance opening")

    intake = json.loads(intake_path.read_text())
    if intake.get("schema") != "ttf_genetic_phylogatr_phase1_intake_receipt_v0.2":
        raise RuntimeError("unexpected regenerated Phase-1 intake schema")
    if intake.get("status") != "PHASE1_FROZEN":
        raise RuntimeError("regenerated Phase-1 intake did not freeze")
    if intake.get("phase1_manifest_sha256") != manifest_sha:
        raise RuntimeError("regenerated Phase-1 intake/manifest linkage drift")
    if intake.get("phase1_geometry_csv_sha256") != geometry_sha:
        raise RuntimeError("regenerated Phase-1 intake/geometry linkage drift")
    for key in (
        "sequence_characters_opened",
        "sequence_characters_hashed",
        "pairwise_genetic_distances_opened",
        "empirical_ttf_opened",
    ):
        if intake.get(key) is not False:
            raise RuntimeError(f"regenerated Phase-1 intake firewall open: {key}")

    projection = intake.get("source_projection")
    if not isinstance(projection, dict):
        raise RuntimeError("regenerated Phase-1 intake lacks source_projection")
    if projection.get("schema") != "ttf_genetic_phylogatr_source_projection_receipt_v0.1":
        raise RuntimeError("unexpected regenerated source-projection schema")
    if projection.get("status") != "ANIMALIA_METADATA_PROJECTION_APPLIED":
        raise RuntimeError("regenerated source projection did not apply")
    if projection.get("raw_archive_sha256") != raw_sha:
        raise RuntimeError("regenerated source projection raw-archive linkage drift")
    expected_projection = rule["expected_projection"]
    for key in (
        "original_row_count",
        "retained_animalia_row_count",
        "removed_non_animalia_row_count",
        "projected_genes_sha256",
        "cite_sha256",
    ):
        if projection.get(key) != expected_projection[key]:
            raise RuntimeError(f"regenerated source projection drift: {key}")
    if projection.get("sequence_identity_opened") is not False:
        raise RuntimeError("regenerated source projection indicates identity opening")

    provenance = manifest.get("provenance", {})
    if provenance.get("genes_sha256") != expected_projection["projected_genes_sha256"]:
        raise RuntimeError("regenerated Phase-1 projected genes provenance drift")
    if provenance.get("cite_sha256") != expected_projection["cite_sha256"]:
        raise RuntimeError("regenerated Phase-1 cite provenance drift")

    return {
        "schema": "ttf_genetic_phylogatr_phase1_source_regeneration_receipt_v0.1",
        "status": "EXACT_PHASE1_MANIFEST_REGENERATED_RESPONSE_BLIND",
        "raw_archive_sha256": raw_sha,
        "raw_archive_size_bytes": int(archive.stat().st_size),
        "phase1_manifest_sha256": manifest_sha,
        "phase1_geometry_csv_sha256": geometry_sha,
        "dataset_digest_sha256": manifest["dataset_digest_sha256"],
        "geometry_fingerprint_sha256": manifest["geometry_fingerprint_sha256"],
        "species_count": int(census["final_species"]),
        "train_count": int(split["train_count"]),
        "eval_count": int(split["eval_count"]),
        "selected_panels_count": int(len(manifest.get("selected_panels", {}))),
        "projected_genes_sha256": projection["projected_genes_sha256"],
        "cite_sha256": projection["cite_sha256"],
        "outcome_firewall": {
            "confirmatory_sequence_identity_opened": False,
            "confirmatory_pairwise_genetic_distances_opened": False,
            "confirmatory_ttf_statistic_opened": False,
            "decker_empirical_genetic_outcomes_opened": False,
        },
        "claim_boundary": (
            "Exact response-blind Phase-1 bytes and selected_panels metadata were "
            "regenerated from the exact frozen raw archive. No nucleotide identity, "
            "pairwise genetic distance, or empirical TTF statistic was opened."
        ),
    }


def main() -> int:
    ap = argparse.ArgumentParser(
        description=(
            "Regenerate the exact lost response-blind Phase-1 manifest from the frozen "
            "authenticated phylogatR archive and verify byte identity before Phase 4."
        )
    )
    ap.add_argument("--archive", type=Path, required=True)
    ap.add_argument("--output-dir", type=Path, required=True)
    ap.add_argument("--repo-root", type=Path, default=Path("."))
    args = ap.parse_args()

    repo_root = args.repo_root.resolve()
    archive = args.archive.resolve()
    output_dir = args.output_dir.resolve()
    rule_path = repo_root / REGENERATION_RULE
    projected_phase1 = repo_root / PROJECTED_PHASE1
    if not archive.is_file():
        raise FileNotFoundError(f"phylogatR archive not found: {archive}")
    if not rule_path.is_file() or not projected_phase1.is_file():
        raise FileNotFoundError("required frozen Phase-1 regeneration files are missing")
    rule = _load_rule(rule_path)

    raw = rule["raw_source"]
    if archive.stat().st_size != int(raw["archive_size_bytes"]):
        raise RuntimeError("raw phylogatR archive size differs from frozen source")
    if sha256_path(archive) != str(raw["raw_archive_sha256"]):
        raise RuntimeError("raw phylogatR archive SHA256 differs from frozen source")

    output_dir.mkdir(parents=True, exist_ok=True)
    receipt_path = output_dir / "phase1_source_regeneration_receipt.json"
    if receipt_path.exists():
        raise FileExistsError(
            "Phase-1 source-regeneration receipt already exists; use a new output directory"
        )

    command = [
        sys.executable,
        str(projected_phase1),
        "--archive",
        str(archive),
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
        raise RuntimeError(
            "response-blind projected Phase-1 regeneration failed.\n"
            + completed.stdout
            + "\n"
            + completed.stderr
        )

    receipt = validate_regenerated_phase1(archive, output_dir, rule)
    receipt["regeneration_rule_sha256"] = sha256_path(rule_path)
    receipt["projected_phase1_entrypoint_sha256"] = sha256_path(projected_phase1)
    receipt_path.write_text(
        json.dumps(receipt, indent=2, sort_keys=True, allow_nan=False) + "\n"
    )
    print(
        json.dumps(
            {
                "status": receipt["status"],
                "phase1_manifest_sha256": receipt["phase1_manifest_sha256"],
                "species_count": receipt["species_count"],
                "selected_panels_count": receipt["selected_panels_count"],
                "sequence_identity_opened": False,
                "pairwise_genetic_distances_opened": False,
                "empirical_ttf_opened": False,
            },
            sort_keys=True,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
