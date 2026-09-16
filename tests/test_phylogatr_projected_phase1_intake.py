from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import zipfile


GENES_HEADERS = (
    "gene",
    "dir",
    "proportion_retained",
    "num_seqs_unaligned",
    "num_seqs_aligned",
    "kingdom",
    "phylum",
    "class",
    "order",
    "family",
    "genus",
    "species",
    "subspecies",
    "different_genbank_species",
)

OCCURRENCE_HEADERS = (
    "phylogatr_id",
    "accession",
    "source_id",
    "latitude",
    "longitude",
    "basis_of_record",
    "coordinate_uncertainty_in_meters",
    "issue",
    "flag",
)


def _write_tsv(path: Path, headers: tuple[str, ...], rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=headers, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def _build_archive(tmp_path: Path) -> Path:
    root = tmp_path / "source" / "phylogatr-results"
    root.mkdir(parents=True)
    (root / "cite.txt").write_text("synthetic source\n", encoding="utf-8")

    genes: list[dict[str, str]] = []
    for index in range(150):
        genus = f"Testgenus{index:03d}"
        epithet = "alpha"
        species = f"{genus} {epithet}"
        stem = f"{genus}-{epithet}-COI"
        rel = f"Insecta/TestOrder/Testidae/{genus}-{epithet}"
        genes.append(
            {
                "gene": stem,
                "dir": rel,
                "proportion_retained": "1.0",
                "num_seqs_unaligned": "12",
                "num_seqs_aligned": "12",
                "kingdom": "Animalia",
                "phylum": "Arthropoda",
                "class": "Insecta",
                "order": "TestOrder",
                "family": "Testidae",
                "genus": genus,
                "species": species,
                "subspecies": "",
                "different_genbank_species": "",
            }
        )
        directory = root / rel
        directory.mkdir(parents=True)
        headers = [f"T{index:03d}_{j:02d}" for j in range(12)]
        with (directory / f"{stem}.afa").open("w", encoding="utf-8") as handle:
            for header in headers:
                handle.write(f">{header}\nACGT\n")
        occurrence_rows = []
        for j, header in enumerate(headers):
            occurrence_rows.append(
                {
                    "phylogatr_id": header,
                    "accession": header,
                    "source_id": str(index * 100 + j),
                    "latitude": str(-40.0 + j * 2.0 + index * 0.0001),
                    "longitude": str(-120.0 + j * 3.0 + index * 0.0001),
                    "basis_of_record": "PRESERVED_SPECIMEN",
                    "coordinate_uncertainty_in_meters": "0",
                    "issue": "",
                    "flag": "",
                }
            )
        _write_tsv(directory / "occurrences.txt", OCCURRENCE_HEADERS, occurrence_rows)

    # A broader-than-requested source row. Deliberately do not create its
    # per-species files: the projection must remove the row before Phase 1.
    genes.append(
        {
            "gene": "Plantus-beta-COI",
            "dir": "Magnoliopsida/TestOrder/Testidae/Plantus-beta",
            "proportion_retained": "1.0",
            "num_seqs_unaligned": "12",
            "num_seqs_aligned": "12",
            "kingdom": "Plantae",
            "phylum": "Tracheophyta",
            "class": "Magnoliopsida",
            "order": "TestOrder",
            "family": "Testidae",
            "genus": "Plantus",
            "species": "Plantus beta",
            "subspecies": "",
            "different_genbank_species": "",
        }
    )
    _write_tsv(root / "genes.txt", GENES_HEADERS, genes)

    archive = tmp_path / "phylogatr_results.zip"
    with zipfile.ZipFile(archive, "w", compression=zipfile.ZIP_DEFLATED) as handle:
        for path in sorted(root.rglob("*")):
            if path.is_file():
                handle.write(path, path.relative_to(root.parent))
    return archive


def test_projected_archive_phase1_keeps_raw_archive_provenance(tmp_path: Path) -> None:
    archive = _build_archive(tmp_path)
    output_dir = tmp_path / "phase1"
    command = [
        sys.executable,
        "scripts/run_phylogatr_projected_phase1_intake.py",
        "--archive",
        str(archive),
        "--output-dir",
        str(output_dir),
        "--repo-root",
        ".",
    ]
    completed = subprocess.run(command, check=False, text=True, capture_output=True)
    assert completed.returncode == 0, completed.stderr

    receipt = json.loads((output_dir / "phase1_intake_receipt.json").read_text())
    manifest = json.loads((output_dir / "phase1_manifest.json").read_text())
    projection = receipt["source_projection"]

    assert receipt["status"] == "PHASE1_FROZEN"
    assert manifest["status"] == "FROZEN_RESPONSE_BLIND_PHASE1_GEOMETRY"
    assert manifest["census"]["genes_rows"] == 150
    assert manifest["census"]["final_species"] == 150
    assert projection["status"] == "ANIMALIA_METADATA_PROJECTION_APPLIED"
    assert projection["raw_archive_sha256"] == hashlib.sha256(archive.read_bytes()).hexdigest()
    assert projection["original_row_count"] == 151
    assert projection["retained_animalia_row_count"] == 150
    assert projection["removed_non_animalia_row_count"] == 1
    assert projection["kingdom_counts"] == {"Animalia": 150, "Plantae": 1}
    assert projection["sequence_identity_opened"] is False
    assert receipt["sequence_characters_opened"] is False
    assert receipt["pairwise_genetic_distances_opened"] is False
    assert receipt["empirical_ttf_opened"] is False
