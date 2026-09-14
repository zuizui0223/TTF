from __future__ import annotations

import csv
import json
import subprocess
import sys
from pathlib import Path

from ttf.phylogatr_confirmatory import GENES_HEADERS, OCCURRENCE_HEADERS


def _write_tsv(path: Path, headers: tuple[str, ...], rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=headers, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def _fresh_archive(root: Path) -> None:
    root.mkdir()
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


def _run(root: Path, output: Path) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        [
            sys.executable,
            "scripts/run_phylogatr_confirmatory_phase1_intake.py",
            "--root",
            str(root),
            "--output-dir",
            str(output),
            "--repo-root",
            ".",
        ],
        check=False,
        text=True,
        capture_output=True,
    )


def test_phase1_intake_stops_before_phase2_and_writes_immutable_receipt(tmp_path: Path) -> None:
    root = tmp_path / "phylogatr-results"
    output = tmp_path / "freeze"
    _fresh_archive(root)

    completed = _run(root, output)
    assert completed.returncode == 0, completed.stderr

    manifest = json.loads((output / "phase1_manifest.json").read_text())
    receipt = json.loads((output / "phase1_intake_receipt.json").read_text())
    assert manifest["status"] == "NOT_EVALUABLE_PHASE1_PANEL_TOO_SMALL"
    assert receipt["status"] == "NOT_EVALUABLE_PHASE1"
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
