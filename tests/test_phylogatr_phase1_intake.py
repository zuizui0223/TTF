from __future__ import annotations

import csv
import json
import stat
import subprocess
import sys
import tarfile
from pathlib import Path
import zipfile

from ttf.phylogatr_confirmatory import GENES_HEADERS, OCCURRENCE_HEADERS


def _write_tsv(path: Path, headers: tuple[str, ...], rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=headers, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def _fresh_archive(root: Path) -> None:
    root.mkdir(parents=True)
    (root / "cite.txt").write_text("fresh test provenance\n")
    raw_dir = "Reptilia/Squamata/Freshus-alpha"
    raw_gene = "Freshus-alpha-COI"
    _write_tsv(
        root / "genes.txt",
        GENES_HEADERS,
        [
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
                "species": "Freshus alpha",
                "subspecies": "",
                "different_genbank_species": "",
            }
        ],
    )
    directory = root / raw_dir
    directory.mkdir(parents=True)
    fasta = directory / f"{raw_gene}.afa"
    fasta.write_bytes(
        b"".join(
            f">ID{i}\n".encode() + (b"ACGT" * 5) + b"\n"
            for i in range(12)
        )
    )
    occurrences = []
    for i in range(12):
        occurrences.append(
            {
                "phylogatr_id": f"ID{i}",
                "accession": f"ACC{i}",
                "source_id": str(i),
                "latitude": str(10.0 + i * 0.4),
                "longitude": str(20.0 + i * 0.7),
                "basis_of_record": "PRESERVED_SPECIMEN",
                "coordinate_uncertainty_in_meters": "",
                "issue": "",
                "flag": "",
            }
        )
    _write_tsv(directory / "occurrences.txt", OCCURRENCE_HEADERS, occurrences)


def _run(source: Path, output: Path, *, archive: bool = False) -> subprocess.CompletedProcess[str]:
    source_flag = "--archive" if archive else "--root"
    return subprocess.run(
        [
            sys.executable,
            "scripts/run_phylogatr_confirmatory_phase1_intake.py",
            source_flag,
            str(source),
            "--output-dir",
            str(output),
            "--repo-root",
            ".",
        ],
        check=False,
        text=True,
        capture_output=True,
    )


def _zip_tree(root: Path, archive: Path) -> None:
    with zipfile.ZipFile(archive, "w", compression=zipfile.ZIP_DEFLATED) as handle:
        for path in sorted(root.rglob("*")):
            if path.is_file():
                handle.write(path, Path("phylogatr-results") / path.relative_to(root))


def _tar_tree(root: Path, archive: Path) -> None:
    with tarfile.open(archive, "w:gz") as handle:
        handle.add(root, arcname="phylogatr-results", recursive=True)


def test_phase1_intake_stops_before_phase2_and_writes_immutable_receipt(tmp_path: Path) -> None:
    root = tmp_path / "phylogatr-results"
    output = tmp_path / "freeze"
    _fresh_archive(root)

    completed = _run(root, output)
    assert completed.returncode == 0, completed.stderr

    manifest = json.loads((output / "phase1_manifest.json").read_text())
    receipt = json.loads((output / "phase1_intake_receipt.json").read_text())
    assert manifest["status"] == "NOT_EVALUABLE_PHASE1_PANEL_TOO_SMALL"
    assert receipt["schema"] == "ttf_genetic_phylogatr_phase1_intake_receipt_v0.2"
    assert receipt["status"] == "NOT_EVALUABLE_PHASE1"
    assert receipt["source_provenance"]["source_type"] == "directory"
    assert receipt["species"] == 1
    assert receipt["sequence_characters_opened"] is False
    assert receipt["sequence_characters_hashed"] is False
    assert receipt["pairwise_genetic_distances_opened"] is False
    assert receipt["empirical_ttf_opened"] is False
    assert receipt["phase2_automatic_execution"] is False
    assert receipt["next_authorized_action"].startswith("STOP")
    assert not (output / "phase2_manifest.json").exists()
    assert not (output / "phase2_geometry.csv").exists()

    rerun = _run(root, output)
    assert rerun.returncode != 0
    assert "already exist" in rerun.stderr


def test_zip_archive_intake_matches_directory_phase1_and_keeps_identity_closed(tmp_path: Path) -> None:
    root = tmp_path / "source" / "phylogatr-results"
    _fresh_archive(root)
    archive = tmp_path / "fresh.zip"
    _zip_tree(root, archive)

    direct_output = tmp_path / "direct"
    archive_output = tmp_path / "from_zip"
    direct = _run(root, direct_output)
    zipped = _run(archive, archive_output, archive=True)
    assert direct.returncode == 0, direct.stderr
    assert zipped.returncode == 0, zipped.stderr

    direct_manifest = json.loads((direct_output / "phase1_manifest.json").read_text())
    archive_manifest = json.loads((archive_output / "phase1_manifest.json").read_text())
    assert archive_manifest["dataset_digest_sha256"] == direct_manifest["dataset_digest_sha256"]
    assert archive_manifest["geometry_csv_sha256"] == direct_manifest["geometry_csv_sha256"]
    assert archive_manifest["geometry_fingerprint_sha256"] == direct_manifest["geometry_fingerprint_sha256"]

    receipt = json.loads((archive_output / "phase1_intake_receipt.json").read_text())
    provenance = receipt["source_provenance"]
    assert provenance["source_type"] == "archive"
    assert provenance["archive_format"] == "zip"
    assert provenance["archive_name"] == "fresh.zip"
    assert len(provenance["archive_sha256"]) == 64
    assert provenance["temporary_extraction"] is True
    assert provenance["temporary_extraction_deleted_after_phase1"] is True
    assert receipt["sequence_characters_opened"] is False
    assert receipt["sequence_characters_hashed"] is False
    assert sorted(path.name for path in archive_output.iterdir()) == [
        "phase1_geometry.csv",
        "phase1_intake_receipt.json",
        "phase1_manifest.json",
    ]


def test_tar_gz_archive_intake_matches_directory_dataset_digest(tmp_path: Path) -> None:
    root = tmp_path / "source" / "phylogatr-results"
    _fresh_archive(root)
    archive = tmp_path / "fresh.tar.gz"
    _tar_tree(root, archive)

    output = tmp_path / "from_tar"
    completed = _run(archive, output, archive=True)
    assert completed.returncode == 0, completed.stderr
    receipt = json.loads((output / "phase1_intake_receipt.json").read_text())
    assert receipt["source_provenance"]["archive_format"] == "tar.gz"
    assert receipt["source_provenance"]["archive_root_relative"] == "phylogatr-results"
    assert receipt["sequence_characters_opened"] is False
    assert receipt["phase2_automatic_execution"] is False


def test_archive_intake_rejects_path_traversal_before_extraction(tmp_path: Path) -> None:
    archive = tmp_path / "unsafe.zip"
    with zipfile.ZipFile(archive, "w") as handle:
        handle.writestr("../escaped.txt", "forbidden")
    output = tmp_path / "freeze"
    completed = _run(archive, output, archive=True)
    assert completed.returncode != 0
    assert "unsafe archive member path" in completed.stderr
    assert not (tmp_path / "escaped.txt").exists()
    assert not (output / "phase1_manifest.json").exists()


def test_archive_intake_rejects_zip_symlink(tmp_path: Path) -> None:
    archive = tmp_path / "symlink.zip"
    with zipfile.ZipFile(archive, "w") as handle:
        info = zipfile.ZipInfo("phylogatr-results/link")
        info.create_system = 3
        info.external_attr = (stat.S_IFLNK | 0o777) << 16
        handle.writestr(info, "../target")
    output = tmp_path / "freeze"
    completed = _run(archive, output, archive=True)
    assert completed.returncode != 0
    assert "symlink archive member is forbidden" in completed.stderr
    assert not (output / "phase1_manifest.json").exists()
