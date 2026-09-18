from __future__ import annotations

import hashlib
import json
import runpy
from pathlib import Path

import pytest


MODULE = runpy.run_path("scripts/regenerate_phylogatr_phase1_for_phase4.py")
validate_regenerated_phase1 = MODULE["validate_regenerated_phase1"]
RULE = Path("docs/supporting/genetic_phylogatr_phase1_source_regeneration_v0.1.json")


def _sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _fixture(tmp_path: Path) -> tuple[Path, Path, dict]:
    archive = tmp_path / "source.zip"
    archive.write_bytes(b"frozen-archive-bytes")
    output = tmp_path / "phase1"
    output.mkdir()

    geometry = output / "phase1_geometry.csv"
    geometry.write_text("species,locality_index\nalpha,0\n")

    projected_genes_sha = "1" * 64
    cite_sha = "2" * 64
    manifest = {
        "schema": "ttf_genetic_phylogatr_confirmatory_phase1_geometry_v0.1",
        "status": "FROZEN_RESPONSE_BLIND_PHASE1_GEOMETRY",
        "dataset_digest_sha256": "3" * 64,
        "geometry_fingerprint_sha256": "4" * 64,
        "census": {"final_species": 1},
        "split": {"train_count": 1, "eval_count": 0},
        "response_blind": {
            "sequence_characters_used": False,
            "sequence_characters_hashed": False,
            "pairwise_genetic_distances_opened": False,
            "ttf_statistic_opened": False,
            "decker_empirical_genetic_outcomes_opened": False,
        },
        "provenance": {
            "genes_sha256": projected_genes_sha,
            "cite_sha256": cite_sha,
        },
        "selected_panels": {"alpha": {"aligned_fasta_relative": "alpha.fasta"}},
        "confirmatory_sequence_identity_opened": False,
        "confirmatory_pairwise_genetic_distances_opened": False,
    }
    manifest_path = output / "phase1_manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")

    projection = {
        "schema": "ttf_genetic_phylogatr_source_projection_receipt_v0.1",
        "status": "ANIMALIA_METADATA_PROJECTION_APPLIED",
        "raw_archive_sha256": _sha(archive),
        "original_row_count": 2,
        "retained_animalia_row_count": 1,
        "removed_non_animalia_row_count": 1,
        "projected_genes_sha256": projected_genes_sha,
        "cite_sha256": cite_sha,
        "sequence_identity_opened": False,
    }
    intake = {
        "schema": "ttf_genetic_phylogatr_phase1_intake_receipt_v0.2",
        "status": "PHASE1_FROZEN",
        "phase1_manifest_sha256": _sha(manifest_path),
        "phase1_geometry_csv_sha256": _sha(geometry),
        "sequence_characters_opened": False,
        "sequence_characters_hashed": False,
        "pairwise_genetic_distances_opened": False,
        "empirical_ttf_opened": False,
        "source_projection": projection,
    }
    (output / "phase1_intake_receipt.json").write_text(
        json.dumps(intake, indent=2, sort_keys=True) + "\n"
    )

    rule = {
        "raw_source": {
            "archive_size_bytes": archive.stat().st_size,
            "raw_archive_sha256": _sha(archive),
        },
        "expected_phase1": {
            "phase1_manifest_sha256": _sha(manifest_path),
            "phase1_geometry_csv_sha256": _sha(geometry),
            "dataset_digest_sha256": manifest["dataset_digest_sha256"],
            "geometry_fingerprint_sha256": manifest["geometry_fingerprint_sha256"],
            "species_count": 1,
            "train_count": 1,
            "eval_count": 0,
        },
        "expected_projection": {
            key: projection[key]
            for key in (
                "original_row_count",
                "retained_animalia_row_count",
                "removed_non_animalia_row_count",
                "projected_genes_sha256",
                "cite_sha256",
            )
        },
    }
    return archive, output, rule


def test_exact_response_blind_phase1_regeneration_validation(tmp_path: Path) -> None:
    archive, output, rule = _fixture(tmp_path)
    receipt = validate_regenerated_phase1(archive, output, rule)
    assert receipt["status"] == "EXACT_PHASE1_MANIFEST_REGENERATED_RESPONSE_BLIND"
    assert receipt["species_count"] == 1
    assert receipt["selected_panels_count"] == 1
    assert all(value is False for value in receipt["outcome_firewall"].values())


def test_phase1_regeneration_rejects_manifest_byte_drift(tmp_path: Path) -> None:
    archive, output, rule = _fixture(tmp_path)
    with (output / "phase1_manifest.json").open("a") as handle:
        handle.write(" ")
    with pytest.raises(RuntimeError, match="not byte-identical"):
        validate_regenerated_phase1(archive, output, rule)


def test_production_regeneration_contract_is_closed_and_pinned() -> None:
    rule = json.loads(RULE.read_text())
    assert rule["status"] == "FROZEN_BEFORE_PHASE4_IDENTITY_OPENING"
    assert rule["raw_source"]["archive_size_bytes"] == 274988692
    assert (
        rule["raw_source"]["raw_archive_sha256"]
        == "5a0fd9ac25893c749d14186fbcce4a46b99163c9d810b36e40eebce7bece61a5"
    )
    assert (
        rule["expected_phase1"]["phase1_manifest_sha256"]
        == "daeedabbeb99fa304568f0f2b2e1c46e0dd6221fbac02b0491faf2435ef7b7b3"
    )
    assert rule["expected_phase1"]["species_count"] == 250
    assert rule["expected_projection"]["retained_animalia_row_count"] == 80495
    assert all(value is False for value in rule["authority"].values())
    assert all(value is False for value in rule["outcome_firewall"].values())

    authorizer = Path("scripts/authorize_phylogatr_phase4_identity_opening.py").read_text()
    for required in (
        "scripts/regenerate_phylogatr_phase1_for_phase4.py",
        "scripts/run_phylogatr_projected_phase1_intake.py",
        "scripts/run_phylogatr_confirmatory_phase1_intake.py",
        "scripts/freeze_phylogatr_confirmatory_phase1.py",
        "docs/supporting/genetic_phylogatr_phase1_source_regeneration_v0.1.json",
    ):
        assert required in authorizer
