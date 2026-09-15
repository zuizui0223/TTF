from __future__ import annotations

import csv
import json
import subprocess
import sys
from pathlib import Path
import zipfile

from ttf.phylogatr_confirmatory import GENES_HEADERS, OCCURRENCE_HEADERS, sha256_path


def _write_tsv(path: Path, headers: tuple[str, ...], rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=headers, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def _fresh_panel(root: Path, *, n_species: int = 150) -> None:
    root.mkdir(parents=True)
    (root / "cite.txt").write_text("fresh Phase-2 archive test provenance\n")
    gene_rows: list[dict[str, str]] = []
    for species_index in range(n_species):
        epithet = f"sp{species_index:03d}"
        species = f"Freshus {epithet}"
        raw_dir = f"Reptilia/Squamata/Freshus-{epithet}"
        raw_gene = f"Freshus-{epithet}-COI"
        gene_rows.append(
            {
                "gene": raw_gene,
                "dir": raw_dir,
                "proportion_retained": "1.0",
                "num_seqs_unaligned": "12",
                "num_seqs_aligned": "12",
                "kingdom": "Animalia",
                "phylum": "Chordata",
                "class": "Reptilia",
                "order": "Squamata",
                "family": "Exampleidae",
                "genus": "Freshus",
                "species": species,
                "subspecies": "",
                "different_genbank_species": "",
            }
        )
        directory = root / raw_dir
        directory.mkdir(parents=True)
        fasta = directory / f"{raw_gene}.afa"
        fasta.write_bytes(
            b"".join(
                f">S{species_index:03d}_ID{i}\n".encode() + (b"ACGT" * 5) + b"\n"
                for i in range(12)
            )
        )
        occurrences = [
            {
                "phylogatr_id": f"S{species_index:03d}_ID{i}",
                "accession": f"ACC{species_index:03d}_{i}",
                "source_id": f"{species_index}_{i}",
                "latitude": str(-30.0 + species_index * 0.01 + i * 0.08),
                "longitude": str(20.0 + i * 0.17),
                "basis_of_record": "PRESERVED_SPECIMEN",
                "coordinate_uncertainty_in_meters": "",
                "issue": "",
                "flag": "",
            }
            for i in range(12)
        ]
        _write_tsv(directory / "occurrences.txt", OCCURRENCE_HEADERS, occurrences)
    _write_tsv(root / "genes.txt", GENES_HEADERS, gene_rows)


def _zip_tree(root: Path, archive: Path) -> None:
    with zipfile.ZipFile(archive, "w", compression=zipfile.ZIP_DEFLATED) as handle:
        for path in sorted(root.rglob("*")):
            if path.is_file():
                handle.write(path, Path("phylogatr-results") / path.relative_to(root))


def _run_phase1(archive: Path, output: Path) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        [
            sys.executable,
            "scripts/run_phylogatr_confirmatory_phase1_intake.py",
            "--archive",
            str(archive),
            "--output-dir",
            str(output),
            "--repo-root",
            ".",
        ],
        check=False,
        text=True,
        capture_output=True,
    )


def _run_phase2(archive: Path, phase1: Path, output: Path) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        [
            sys.executable,
            "scripts/run_phylogatr_confirmatory_phase2_intake.py",
            "--archive",
            str(archive),
            "--phase1-dir",
            str(phase1),
            "--output-dir",
            str(output),
            "--repo-root",
            ".",
        ],
        check=False,
        text=True,
        capture_output=True,
    )


def _fake_phase1_chain(phase1_dir: Path, archive: Path) -> None:
    phase1_dir.mkdir(parents=True)
    geometry = phase1_dir / "phase1_geometry.csv"
    geometry.write_text("species,locality_index,x_km,y_km,z_km,graph_k\n")
    manifest = {
        "schema": "ttf_genetic_phylogatr_confirmatory_phase1_geometry_v0.1",
        "status": "FROZEN_RESPONSE_BLIND_PHASE1_GEOMETRY",
        "dataset_digest_sha256": "a" * 64,
        "geometry_fingerprint_sha256": "b" * 64,
        "confirmatory_sequence_identity_opened": False,
        "confirmatory_pairwise_genetic_distances_opened": False,
    }
    manifest_path = phase1_dir / "phase1_manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    receipt = {
        "schema": "ttf_genetic_phylogatr_phase1_intake_receipt_v0.2",
        "status": "PHASE1_FROZEN",
        "source_provenance": {
            "source_type": "archive",
            "archive_sha256": sha256_path(archive),
            "archive_format": "zip",
            "archive_root_relative": "phylogatr-results",
            "cite_sha256": "c" * 64,
            "genes_sha256": "d" * 64,
        },
        "phase1_manifest_sha256": sha256_path(manifest_path),
        "phase1_geometry_csv_sha256": sha256_path(geometry),
        "dataset_digest_sha256": "a" * 64,
        "geometry_fingerprint_sha256": "b" * 64,
    }
    (phase1_dir / "phase1_intake_receipt.json").write_text(
        json.dumps(receipt, indent=2, sort_keys=True) + "\n"
    )


def test_archive_direct_phase2_reuses_exact_phase1_archive_and_keeps_identity_closed(
    tmp_path: Path,
) -> None:
    root = tmp_path / "source" / "phylogatr-results"
    _fresh_panel(root)
    archive = tmp_path / "fresh.zip"
    _zip_tree(root, archive)

    phase1 = tmp_path / "phase1"
    phase2 = tmp_path / "phase2"
    completed1 = _run_phase1(archive, phase1)
    assert completed1.returncode == 0, completed1.stderr
    receipt1 = json.loads((phase1 / "phase1_intake_receipt.json").read_text())
    assert receipt1["status"] == "PHASE1_FROZEN"
    assert receipt1["species"] == 150

    completed2 = _run_phase2(archive, phase1, phase2)
    assert completed2.returncode == 0, completed2.stderr
    manifest2 = json.loads((phase2 / "phase2_manifest.json").read_text())
    receipt2 = json.loads((phase2 / "phase2_intake_receipt.json").read_text())
    assert manifest2["status"] == "PASS_TO_SYNTHETIC_GATE"
    assert manifest2["species"]["survivors"] == 150
    assert receipt2["status"] == "PHASE2_FROZEN_TO_SYNTHETIC_GATE"
    assert receipt2["source_provenance"]["source_type"] == "archive"
    assert receipt2["source_provenance"]["archive_sha256"] == receipt1["source_provenance"]["archive_sha256"]
    assert receipt2["source_provenance"]["temporary_extraction_deleted_after_phase2"] is True
    assert receipt2["character_mask_opened"] is True
    assert receipt2["sequence_identity_opened"] is False
    assert receipt2["pairwise_genetic_distances_opened"] is False
    assert receipt2["empirical_ttf_opened"] is False
    assert receipt2["phase3_automatic_execution"] is False
    assert sorted(path.name for path in phase2.iterdir()) == [
        "phase2_geometry.csv",
        "phase2_intake_receipt.json",
        "phase2_manifest.json",
    ]

    rerun = _run_phase2(archive, phase1, phase2)
    assert rerun.returncode != 0
    assert "already exist" in rerun.stderr


def test_phase2_rejects_archive_changed_after_phase1(tmp_path: Path) -> None:
    root = tmp_path / "source" / "phylogatr-results"
    _fresh_panel(root)
    archive = tmp_path / "fresh.zip"
    _zip_tree(root, archive)
    phase1 = tmp_path / "phase1"
    completed1 = _run_phase1(archive, phase1)
    assert completed1.returncode == 0, completed1.stderr

    with zipfile.ZipFile(archive, "a") as handle:
        handle.writestr("extra_after_phase1.txt", "drift")
    phase2 = tmp_path / "phase2"
    completed2 = _run_phase2(archive, phase1, phase2)
    assert completed2.returncode != 0
    assert "archive SHA256 differs" in completed2.stderr
    assert not (phase2 / "phase2_manifest.json").exists()


def test_phase2_archive_intake_rejects_path_traversal_before_extraction(tmp_path: Path) -> None:
    archive = tmp_path / "unsafe.zip"
    with zipfile.ZipFile(archive, "w") as handle:
        handle.writestr("../escaped.txt", "forbidden")
    phase1 = tmp_path / "phase1"
    _fake_phase1_chain(phase1, archive)
    phase2 = tmp_path / "phase2"
    completed = _run_phase2(archive, phase1, phase2)
    assert completed.returncode != 0
    assert "unsafe archive member path" in completed.stderr
    assert not (tmp_path / "escaped.txt").exists()
    assert not (phase2 / "phase2_manifest.json").exists()
