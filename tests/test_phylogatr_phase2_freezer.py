from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path
import subprocess
import sys

import numpy as np

from ttf.genetic_geometry import prepare_genetic_sampling_geometry
from ttf.phylogatr_confirmatory import (
    OCCURRENCE_HEADERS,
    canonical_latlon_for_geometry,
    fasta_header_sha256,
    fasta_headers_only,
    latlon_to_ecef_km,
    sha256_path,
)


def _write_occurrences(path: Path, latlon: np.ndarray) -> None:
    rows = []
    for index, (lat, lon) in enumerate(latlon):
        rows.append(
            {
                "phylogatr_id": f"H{index}",
                "accession": f"ACC{index}",
                "source_id": str(index),
                "latitude": str(float(lat)),
                "longitude": str(float(lon)),
                "basis_of_record": "PRESERVED_SPECIMEN",
                "coordinate_uncertainty_in_meters": "",
                "issue": "",
                "flag": "",
            }
        )
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=OCCURRENCE_HEADERS, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def _write_fasta(path: Path, sequences: list[str]) -> None:
    lines: list[str] = []
    for index, sequence in enumerate(sequences):
        lines.extend([f">H{index}\n", sequence + "\n"])
    path.write_text("".join(lines))


def _make_species(root: Path, species: str, *, fail_mask: bool) -> tuple[dict, list[dict[str, str]]]:
    stem = species.replace(" ", "-")
    directory = root / "Reptilia" / "Squamata" / stem
    directory.mkdir(parents=True)
    fasta = directory / f"{stem}-COI.afa"
    occurrence = directory / "occurrences.txt"
    source_latlon = np.asarray(
        [[0.0, 10.0], [1.0, 11.0], [2.0, 12.0], [3.0, 13.0]], dtype=float
    )
    _write_occurrences(occurrence, source_latlon)
    sequences = ["ACGTACGTAC"] * 4
    if fail_mask:
        sequences[-1] = "NNNNNNNNNN"
    _write_fasta(fasta, sequences)

    raw_ecef = latlon_to_ecef_km(source_latlon)
    geometry = prepare_genetic_sampling_geometry(raw_ecef, k=3)
    canonical_latlon = canonical_latlon_for_geometry(source_latlon, geometry)
    csv_rows: list[dict[str, str]] = []
    for index, ((lat, lon), (x, y, z)) in enumerate(
        zip(canonical_latlon, geometry.coordinates)
    ):
        csv_rows.append(
            {
                "species": species,
                "locality_index": str(index),
                "latitude": repr(float(lat)),
                "longitude": repr(float(lon)),
                "x_km": repr(float(x)),
                "y_km": repr(float(y)),
                "z_km": repr(float(z)),
                "graph_k": "3",
            }
        )
    meta = {
        "aligned_fasta_relative": str(fasta.relative_to(root)),
        "occurrence_relative": str(occurrence.relative_to(root)),
        "aligned_header_sha256": fasta_header_sha256(fasta_headers_only(fasta)),
        "occurrence_sha256": sha256_path(occurrence),
        "unique_localities": 4,
        "graph_k": 3,
        "edge_count": int(geometry.n_edges),
    }
    return meta, csv_rows


def _sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _build_phase1_fixture(tmp_path: Path) -> tuple[Path, Path, Path]:
    root = tmp_path / "phylogatr-results"
    root.mkdir()
    (root / "cite.txt").write_text("fresh-fixture\n")
    (root / "genes.txt").write_text("frozen-fixture-schema-placeholder\n")

    survivor_meta, survivor_rows = _make_species(root, "Freshus alpha", fail_mask=False)
    failure_meta, failure_rows = _make_species(root, "Freshus beta", fail_mask=True)
    csv_path = tmp_path / "phase1.csv"
    fieldnames = [
        "species",
        "locality_index",
        "latitude",
        "longitude",
        "x_km",
        "y_km",
        "z_km",
        "graph_k",
    ]
    with csv_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(survivor_rows + failure_rows)

    manifest = {
        "schema": "ttf_genetic_phylogatr_confirmatory_phase1_geometry_v0.1",
        "status": "FROZEN_RESPONSE_BLIND_PHASE1_GEOMETRY",
        "dataset_digest_sha256": "fixture-dataset-digest",
        "geometry_fingerprint_sha256": "fixture-phase1-fingerprint",
        "geometry_csv_sha256": _sha(csv_path),
        "provenance": {
            "cite_sha256": _sha(root / "cite.txt"),
            "genes_sha256": _sha(root / "genes.txt"),
        },
        "response_blind": {
            "sequence_characters_used": False,
            "sequence_characters_hashed": False,
        },
        "confirmatory_sequence_identity_opened": False,
        "confirmatory_pairwise_genetic_distances_opened": False,
        "census": {"final_species": 2},
        "split": {
            "train_species": ["Freshus alpha"],
            "eval_species": ["Freshus beta"],
            "train_count": 1,
            "eval_count": 1,
        },
        "selected_panels": {
            "Freshus alpha": survivor_meta,
            "Freshus beta": failure_meta,
        },
    }
    manifest_path = tmp_path / "phase1.json"
    manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    return root, csv_path, manifest_path


def test_phase2_freezer_filters_species_without_rewire_or_resplit(tmp_path: Path) -> None:
    root, phase1_csv, phase1_manifest = _build_phase1_fixture(tmp_path)
    out_csv = tmp_path / "phase2.csv"
    out_json = tmp_path / "phase2.json"
    command = [
        sys.executable,
        "scripts/freeze_phylogatr_confirmatory_phase2_mask.py",
        "--root",
        str(root),
        "--phase1-manifest",
        str(phase1_manifest),
        "--phase1-csv",
        str(phase1_csv),
        "--phase2-rule",
        "docs/supporting/genetic_phylogatr_phase2_mask_rule_v0.1.json",
        "--output-csv",
        str(out_csv),
        "--output-manifest",
        str(out_json),
    ]
    completed = subprocess.run(command, check=False, text=True, capture_output=True)
    assert completed.returncode == 0, completed.stderr
    payload = json.loads(out_json.read_text())
    assert payload["status"] == "NOT_EVALUABLE_PHASE2_PANEL_TOO_SMALL"
    assert payload["species"] == {
        "phase1": 2,
        "survivors": 1,
        "failed_character_support": 1,
        "minimum_required": 150,
    }
    assert payload["species_ledger"]["Freshus alpha"]["survives"] is True
    assert payload["species_ledger"]["Freshus beta"]["survives"] is False
    assert payload["species_ledger"]["Freshus beta"]["invalid_edges"] > 0
    assert payload["split"]["inherit_phase1_without_resplitting"] is True
    assert payload["split"]["train_species"] == ["Freshus alpha"]
    assert payload["split"]["eval_species"] == []
    assert payload["character_mask_opened"] is True
    assert payload["confirmatory_sequence_identity_opened"] is False
    assert payload["confirmatory_pairwise_genetic_distances_opened"] is False
    assert payload["confirmatory_ttf_statistic_opened"] is False

    with out_csv.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    assert rows
    assert {row["species"] for row in rows} == {"Freshus alpha"}
    # Graph is not rewired: all four frozen localities survive unchanged.
    assert [int(row["locality_index"]) for row in rows] == [0, 1, 2, 3]


def test_phase2_source_hash_drift_is_global_stop(tmp_path: Path) -> None:
    root, phase1_csv, phase1_manifest = _build_phase1_fixture(tmp_path)
    payload = json.loads(phase1_manifest.read_text())
    occurrence = root / payload["selected_panels"]["Freshus alpha"]["occurrence_relative"]
    occurrence.write_text(occurrence.read_text() + "\n")
    completed = subprocess.run(
        [
            sys.executable,
            "scripts/freeze_phylogatr_confirmatory_phase2_mask.py",
            "--root",
            str(root),
            "--phase1-manifest",
            str(phase1_manifest),
            "--phase1-csv",
            str(phase1_csv),
            "--phase2-rule",
            "docs/supporting/genetic_phylogatr_phase2_mask_rule_v0.1.json",
            "--output-csv",
            str(tmp_path / "out.csv"),
            "--output-manifest",
            str(tmp_path / "out.json"),
        ],
        check=False,
        text=True,
        capture_output=True,
    )
    assert completed.returncode != 0
    assert "occurrence source drift since phase-1 freeze" in completed.stderr
    assert not (tmp_path / "out.json").exists()
