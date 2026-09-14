from __future__ import annotations

import csv
import json
import subprocess
import sys
from pathlib import Path

import pytest

from ttf.genetic_geometry_io import load_frozen_genetic_geometry_csv, sha256_path
from ttf.geometry import SpeciesGeometry, geometry_fingerprint
from ttf.phylogatr_phase3 import derive_phase3_master_seed, load_phylogatr_phase3_context


FIELDS = ["species", "locality_index", "x_km", "y_km", "z_km", "graph_k"]
RULE = Path("docs/supporting/genetic_phylogatr_phase3_gate_d_rule_v0.1.json")


def _build_phase2_fixture(
    tmp_path: Path,
) -> tuple[Path, Path, Path, tuple[str, ...], tuple[str, ...]]:
    geometry_csv = tmp_path / "phase2_geometry.csv"
    species_names = tuple(f"species_{index:03d}" for index in range(150))
    rows: list[dict[str, object]] = []
    for species_index, species in enumerate(species_names):
        for locality in range(12):
            rows.append(
                {
                    "species": species,
                    "locality_index": locality,
                    "x_km": species_index * 1000.0 + locality * 10.0,
                    "y_km": float(locality % 3),
                    "z_km": float(locality % 5),
                    "graph_k": 2,
                }
            )
    with geometry_csv.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDS)
        writer.writeheader()
        writer.writerows(rows)

    table = load_frozen_genetic_geometry_csv(geometry_csv)
    fingerprint = geometry_fingerprint(
        [
            SpeciesGeometry(species=name, coordinates=table.geometries[name].coordinates)
            for name in table.species
        ]
    )
    geometry_sha = sha256_path(geometry_csv)
    train = species_names[:75]
    evaluation = species_names[75:]
    dataset_digest = "a" * 64

    phase1 = {
        "schema": "ttf_genetic_phylogatr_confirmatory_phase1_geometry_v0.1",
        "status": "FROZEN_RESPONSE_BLIND_PHASE1_GEOMETRY",
        "dataset_digest_sha256": dataset_digest,
        "geometry_csv_sha256": geometry_sha,
        "geometry_fingerprint_sha256": fingerprint,
        "response_blind": {
            "sequence_characters_used": False,
            "sequence_characters_hashed": False,
            "pairwise_genetic_distances_opened": False,
            "ttf_statistic_opened": False,
            "decker_empirical_genetic_outcomes_opened": False,
            "allowed_fasta_information": "header strings only",
        },
        "census": {
            "final_species": 150,
            "minimum_panel_pass": True,
        },
        "split": {
            "train_species": list(train),
            "eval_species": list(evaluation),
            "train_count": 75,
            "eval_count": 75,
        },
        "confirmatory_sequence_identity_opened": False,
        "confirmatory_pairwise_genetic_distances_opened": False,
        "qualification_claim_made": False,
    }
    phase1_path = tmp_path / "phase1.json"
    phase1_path.write_text(json.dumps(phase1, indent=2, sort_keys=True) + "\n")

    phase2 = {
        "schema": "ttf_genetic_phylogatr_confirmatory_phase2_mask_v0.1",
        "status": "PASS_TO_SYNTHETIC_GATE",
        "phase1": {
            "dataset_digest_sha256": dataset_digest,
            "geometry_fingerprint_sha256": fingerprint,
            "geometry_csv_sha256": geometry_sha,
            "species_count": 150,
            "train_count": 75,
            "eval_count": 75,
        },
        "rule_sha256": "d" * 64,
        "mask_contract": {
            "canonical_valid_characters": ["A", "C", "G", "T", "a", "c", "g", "t"],
            "minimum_comparable_fraction": 0.5,
            "nucleotide_identity_persisted": False,
            "pairwise_nucleotide_differences_computed": False,
        },
        "species": {
            "phase1": 150,
            "survivors": 150,
            "failed_character_support": 0,
            "minimum_required": 150,
        },
        "split": {
            "inherit_phase1_without_resplitting": True,
            "train_species": list(train),
            "eval_species": list(evaluation),
            "train_count": 75,
            "eval_count": 75,
        },
        "geometry_csv_sha256": geometry_sha,
        "geometry_fingerprint_sha256": fingerprint,
        "geometry_fingerprint_status": "computed_core_fingerprint",
        "species_ledger": {},
        "character_mask_opened": True,
        "confirmatory_sequence_identity_opened": False,
        "confirmatory_pairwise_genetic_distances_opened": False,
        "confirmatory_ttf_statistic_opened": False,
        "qualification_claim_made": False,
        "next_gate": "repeat_full_dataset_specific_profiled_private_Gate_D_on_phase2_survivor_geometry",
    }
    phase2_path = tmp_path / "phase2.json"
    phase2_path.write_text(json.dumps(phase2, indent=2, sort_keys=True) + "\n")
    return geometry_csv, phase1_path, phase2_path, train, evaluation


def _authorize(tmp_path: Path) -> tuple[Path, Path, Path, Path]:
    geometry_csv, phase1_path, phase2_path, _, _ = _build_phase2_fixture(tmp_path)
    authorization = tmp_path / "phase3_authorization.json"
    command = [
        sys.executable,
        "scripts/authorize_phylogatr_phase3_gate_d.py",
        "--geometry",
        str(geometry_csv),
        "--phase1-manifest",
        str(phase1_path),
        "--phase2-manifest",
        str(phase2_path),
        "--phase3-rule",
        str(RULE),
        "--repo-root",
        ".",
        "--output",
        str(authorization),
    ]
    completed = subprocess.run(command, check=False, text=True, capture_output=True)
    assert completed.returncode == 0, completed.stderr
    return geometry_csv, phase1_path, phase2_path, authorization


def test_phase3_authorization_and_runtime_contract_pass(tmp_path: Path) -> None:
    geometry_csv, phase1_path, phase2_path, authorization = _authorize(tmp_path)
    payload = json.loads(authorization.read_text())
    phase2 = json.loads(phase2_path.read_text())
    expected_seed = derive_phase3_master_seed(
        phase2["phase1"]["dataset_digest_sha256"],
        phase2["geometry_fingerprint_sha256"],
    )
    assert payload["status"] == "authorize_frozen_fresh_phylogatr_phase3_gate_d"
    assert payload["phase1_manifest_sha256"] == sha256_path(phase1_path)
    assert payload["phase2_manifest_sha256"] == sha256_path(phase2_path)
    assert payload["species"]["survivors"] == 150
    assert payload["species"]["train"] == 75
    assert payload["species"]["eval"] == 75
    assert payload["master_seed"] == expected_seed
    assert payload["execution"]["reference_expected_jobs"] == 64
    assert payload["execution"]["observed_expected_jobs"] == 10
    assert all(value is False for value in payload["outcome_firewall"].values())

    context = load_phylogatr_phase3_context(
        geometry_csv,
        phase2_path,
        RULE,
        authorization,
        verify_code=True,
    )
    assert len(context.geometries) == 150
    assert len(context.train_species) == 75
    assert len(context.eval_species) == 75
    assert context.master_seed == expected_seed


def test_phase3_shard_plan_exactly_covers_frozen_replicates(tmp_path: Path) -> None:
    _, _, _, authorization = _authorize(tmp_path)
    output = tmp_path / "plan.json"
    completed = subprocess.run(
        [
            sys.executable,
            "scripts/plan_phylogatr_phase3_shards.py",
            "--phase3-rule",
            str(RULE),
            "--authorization",
            str(authorization),
            "--output",
            str(output),
        ],
        check=False,
        text=True,
        capture_output=True,
    )
    assert completed.returncode == 0, completed.stderr
    plan = json.loads(output.read_text())
    assert plan["reference_job_count"] == 64
    assert plan["observed_job_count"] == 10

    by_config: dict[str, list[int]] = {}
    for job in plan["reference_jobs"]:
        by_config.setdefault(job["configuration"], []).extend(range(job["start"], job["stop"]))
    assert len(by_config) == 8
    assert all(indices == list(range(1999)) for indices in by_config.values())

    by_cell: dict[tuple[float, float], list[int]] = {}
    for job in plan["observed_jobs"]:
        key = (job["shared_fraction"], job["residual_amplitude"])
        by_cell.setdefault(key, []).extend(range(job["start"], job["stop"]))
    assert len(by_cell) == 5
    assert all(indices == list(range(500)) for indices in by_cell.values())


def test_phase3_runtime_rejects_geometry_drift_after_authorization(tmp_path: Path) -> None:
    geometry_csv, _, phase2_path, authorization = _authorize(tmp_path)
    geometry_csv.write_text(geometry_csv.read_text() + "\n")
    with pytest.raises(RuntimeError, match="geometry CSV SHA256 drift"):
        load_phylogatr_phase3_context(
            geometry_csv,
            phase2_path,
            RULE,
            authorization,
            verify_code=False,
        )


def test_phase3_authorization_rejects_identity_opening(tmp_path: Path) -> None:
    geometry_csv, phase1_path, phase2_path, _, _ = _build_phase2_fixture(tmp_path)
    payload = json.loads(phase2_path.read_text())
    payload["confirmatory_sequence_identity_opened"] = True
    phase2_path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    authorization = tmp_path / "forbidden.json"
    command = [
        sys.executable,
        "scripts/authorize_phylogatr_phase3_gate_d.py",
        "--geometry",
        str(geometry_csv),
        "--phase1-manifest",
        str(phase1_path),
        "--phase2-manifest",
        str(phase2_path),
        "--phase3-rule",
        str(RULE),
        "--repo-root",
        ".",
        "--output",
        str(authorization),
    ]
    completed = subprocess.run(command, check=False, text=True, capture_output=True)
    assert completed.returncode != 0
    assert "nucleotide identity opening" in completed.stderr
    assert not authorization.exists()


def test_phase3_authorization_rejects_phase1_phase2_provenance_drift(tmp_path: Path) -> None:
    geometry_csv, phase1_path, phase2_path, _, _ = _build_phase2_fixture(tmp_path)
    phase1 = json.loads(phase1_path.read_text())
    phase1["dataset_digest_sha256"] = "f" * 64
    phase1_path.write_text(json.dumps(phase1, indent=2, sort_keys=True) + "\n")
    authorization = tmp_path / "drift.json"
    completed = subprocess.run(
        [
            sys.executable,
            "scripts/authorize_phylogatr_phase3_gate_d.py",
            "--geometry",
            str(geometry_csv),
            "--phase1-manifest",
            str(phase1_path),
            "--phase2-manifest",
            str(phase2_path),
            "--phase3-rule",
            str(RULE),
            "--repo-root",
            ".",
            "--output",
            str(authorization),
        ],
        check=False,
        text=True,
        capture_output=True,
    )
    assert completed.returncode != 0
    assert "Phase-1/Phase-2 provenance drift" in completed.stderr
    assert not authorization.exists()
