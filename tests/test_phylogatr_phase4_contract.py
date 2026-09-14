from __future__ import annotations

import csv
import json
import subprocess
import sys
from pathlib import Path

from ttf.genetic_geometry_io import load_frozen_genetic_geometry_csv, sha256_path
from ttf.geometry import SpeciesGeometry, geometry_fingerprint
from ttf.phylogatr_phase4 import load_phylogatr_phase4_context


PHASE3_RULE = Path("docs/supporting/genetic_phylogatr_phase3_gate_d_rule_v0.1.json")
SELF_RULE = Path("docs/supporting/genetic_phylogatr_phase3_self_detectability_rule_v0.1.json")
PHASE4_RULE = Path("docs/supporting/genetic_phylogatr_phase4_response_rule_v0.1.json")
OPENING_STATE = Path("benchmarks/frozen/genetic_empirical_opening_state_v0.1.json")
FIELDS = [
    "species",
    "locality_index",
    "latitude",
    "longitude",
    "x_km",
    "y_km",
    "z_km",
    "graph_k",
]


def _build_chain(tmp_path: Path, *, phase3_pass: bool = True, self_pass: bool = False) -> dict[str, Path]:
    geometry = tmp_path / "phase2_geometry.csv"
    names = tuple(f"species_{index:03d}" for index in range(150))
    with geometry.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDS)
        writer.writeheader()
        for species_index, species in enumerate(names):
            for locality in range(12):
                writer.writerow(
                    {
                        "species": species,
                        "locality_index": locality,
                        "latitude": -40.0 + species_index * 0.01 + locality * 0.001,
                        "longitude": 100.0 + locality * 0.01,
                        "x_km": species_index * 1000.0 + locality * 10.0,
                        "y_km": float(locality % 3),
                        "z_km": float(locality % 5),
                        "graph_k": 2,
                    }
                )
    table = load_frozen_genetic_geometry_csv(geometry)
    fingerprint = geometry_fingerprint(
        [
            SpeciesGeometry(species=name, coordinates=table.geometries[name].coordinates)
            for name in table.species
        ]
    )
    geometry_sha = sha256_path(geometry)
    train, evaluation = names[:75], names[75:]
    dataset_digest = "a" * 64

    phase1_payload = {
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
        },
        "provenance": {"cite_sha256": "1" * 64, "genes_sha256": "2" * 64},
        "census": {"final_species": 150, "minimum_panel_pass": True},
        "split": {
            "train_species": list(train),
            "eval_species": list(evaluation),
            "train_count": 75,
            "eval_count": 75,
        },
        "selected_panels": {name: {} for name in names},
        "confirmatory_sequence_identity_opened": False,
        "confirmatory_pairwise_genetic_distances_opened": False,
    }
    phase1 = tmp_path / "phase1.json"
    phase1.write_text(json.dumps(phase1_payload, indent=2, sort_keys=True) + "\n")

    phase2_payload = {
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
            "nucleotide_identity_persisted": False,
            "pairwise_nucleotide_differences_computed": False,
        },
        "species": {"survivors": 150, "minimum_required": 150},
        "split": {
            "inherit_phase1_without_resplitting": True,
            "train_species": list(train),
            "eval_species": list(evaluation),
            "train_count": 75,
            "eval_count": 75,
        },
        "geometry_csv_sha256": geometry_sha,
        "geometry_fingerprint_sha256": fingerprint,
        "species_ledger": {},
        "character_mask_opened": True,
        "confirmatory_sequence_identity_opened": False,
        "confirmatory_pairwise_genetic_distances_opened": False,
        "confirmatory_ttf_statistic_opened": False,
    }
    phase2 = tmp_path / "phase2.json"
    phase2.write_text(json.dumps(phase2_payload, indent=2, sort_keys=True) + "\n")

    phase3_auth = tmp_path / "phase3_auth.json"
    completed = subprocess.run(
        [
            sys.executable,
            "scripts/authorize_phylogatr_phase3_gate_d.py",
            "--geometry", str(geometry),
            "--phase1-manifest", str(phase1),
            "--phase2-manifest", str(phase2),
            "--phase3-rule", str(PHASE3_RULE),
            "--repo-root", ".",
            "--output", str(phase3_auth),
        ],
        check=False,
        text=True,
        capture_output=True,
    )
    assert completed.returncode == 0, completed.stderr

    references = tmp_path / "references.json"
    references.write_text(
        json.dumps(
            {
                "schema": "ttf_genetic_phylogatr_phase3_references_v0.1",
                "status": "ordered_private_reference_families_complete",
                "geometry_fingerprint_sha256": fingerprint,
                "references": {},
                "confirmatory_sequence_identity_opened": False,
                "confirmatory_pairwise_genetic_distances_opened": False,
            },
            sort_keys=True,
        ) + "\n"
    )
    qualification = tmp_path / "qualification.json"
    qualification.write_text(
        json.dumps(
            {
                "schema": "ttf_genetic_phylogatr_phase3_qualification_v0.1",
                "status": "PASS" if phase3_pass else "NOT_EVALUABLE",
                "geometry_fingerprint_sha256": fingerprint,
                "passed": phase3_pass,
                "phase4_identity_opening_eligible": phase3_pass,
                "type1_gate": {"pass": phase3_pass},
                "power_gate": {"pass": phase3_pass},
                "confirmatory_sequence_identity_opened": False,
                "confirmatory_pairwise_genetic_distances_opened": False,
                "confirmatory_ttf_statistic_opened": False,
            },
            sort_keys=True,
        ) + "\n"
    )
    self_refs = tmp_path / "self_refs.json"
    self_refs.write_text(
        json.dumps(
            {
                "schema": "ttf_genetic_phylogatr_phase3_self_references_v0.1",
                "status": "complete_independent_null_reference",
                "geometry_fingerprint_sha256": fingerprint,
                "statistics": [0.0] * 1000,
                "confirmatory_sequence_identity_opened": False,
                "confirmatory_pairwise_genetic_distances_opened": False,
            },
            sort_keys=True,
        ) + "\n"
    )
    self_qualification = tmp_path / "self_qualification.json"
    self_qualification.write_text(
        json.dumps(
            {
                "schema": "ttf_genetic_phylogatr_phase3_self_qualification_v0.1",
                "status": "PASS" if self_pass else "SELF_DETECTABILITY_NOT_QUALIFIED",
                "geometry_fingerprint_sha256": fingerprint,
                "passed": self_pass,
                "confirmatory_sequence_identity_opened": False,
                "confirmatory_pairwise_genetic_distances_opened": False,
                "confirmatory_ttf_statistic_opened": False,
            },
            sort_keys=True,
        ) + "\n"
    )
    return {
        "geometry": geometry,
        "phase1": phase1,
        "phase2": phase2,
        "phase3_auth": phase3_auth,
        "references": references,
        "qualification": qualification,
        "self_refs": self_refs,
        "self_qualification": self_qualification,
    }


def _authorize_phase4(tmp_path: Path, chain: dict[str, Path]) -> tuple[subprocess.CompletedProcess[str], Path]:
    output = tmp_path / "phase4_auth.json"
    completed = subprocess.run(
        [
            sys.executable,
            "scripts/authorize_phylogatr_phase4_identity_opening.py",
            "--geometry", str(chain["geometry"]),
            "--phase1-manifest", str(chain["phase1"]),
            "--phase2-manifest", str(chain["phase2"]),
            "--phase3-rule", str(PHASE3_RULE),
            "--phase3-authorization", str(chain["phase3_auth"]),
            "--references", str(chain["references"]),
            "--qualification", str(chain["qualification"]),
            "--self-rule", str(SELF_RULE),
            "--self-references", str(chain["self_refs"]),
            "--self-qualification", str(chain["self_qualification"]),
            "--phase4-rule", str(PHASE4_RULE),
            "--opening-state", str(OPENING_STATE),
            "--repo-root", ".",
            "--output", str(output),
        ],
        check=False,
        text=True,
        capture_output=True,
    )
    return completed, output


def test_phase4_authorization_allows_complete_self_failure_without_rescuing_negative(tmp_path: Path) -> None:
    chain = _build_chain(tmp_path, phase3_pass=True, self_pass=False)
    completed, authorization = _authorize_phase4(tmp_path, chain)
    assert completed.returncode == 0, completed.stderr
    payload = json.loads(authorization.read_text())
    assert payload["status"] == "AUTHORIZE_EXACT_FRESH_NUCLEOTIDE_IDENTITY_OPENING"
    assert payload["fresh_self_detectability"]["passed"] is False
    assert payload["fresh_self_detectability"]["negative_interpretation_if_not_passed"] == "NOT_EVALUABLE_FOR_LINEAGE_CONDITIONING"
    assert all(value is False for value in payload["outcome_firewall_before_execution"].values())

    context = load_phylogatr_phase4_context(
        chain["geometry"],
        chain["phase1"],
        chain["phase2"],
        PHASE3_RULE,
        chain["phase3_auth"],
        chain["references"],
        chain["qualification"],
        SELF_RULE,
        chain["self_refs"],
        chain["self_qualification"],
        PHASE4_RULE,
        authorization,
        OPENING_STATE,
        verify_code=True,
    )
    assert context.self_qualification["passed"] is False
    assert len(context.geometries) == 150


def test_phase4_authorization_rejects_failed_phase3_gate(tmp_path: Path) -> None:
    chain = _build_chain(tmp_path, phase3_pass=False, self_pass=True)
    completed, authorization = _authorize_phase4(tmp_path, chain)
    assert completed.returncode != 0
    assert "Phase-3 Gate-D did not PASS" in completed.stderr
    assert not authorization.exists()
