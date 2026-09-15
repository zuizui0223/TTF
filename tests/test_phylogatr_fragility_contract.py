from __future__ import annotations

import csv
import json
from pathlib import Path

import pytest

from ttf.genetic_geometry_io import load_frozen_genetic_geometry_csv, sha256_path
from ttf.geometry import SpeciesGeometry, geometry_fingerprint
from ttf.heterogeneous_inference import density_scaled_k
from ttf.phylogatr_fragility import (
    build_fragility_plan,
    retained_locality_indices,
    write_fragility_plan,
)


PHASE3_RULE = Path("docs/supporting/genetic_phylogatr_phase3_gate_d_rule_v0.1.json")
FRAGILITY_RULE = Path(
    "docs/supporting/genetic_phylogatr_gate_d_fragility_diagnostic_v0.1.json"
)
FIELDS = ["species", "locality_index", "x_km", "y_km", "z_km", "graph_k"]


def _fixture(tmp_path: Path) -> tuple[Path, Path]:
    geometry_csv = tmp_path / "phase2_geometry.csv"
    species_names = tuple(f"fresh_species_{index:03d}" for index in range(150))
    n_localities = 24
    graph_k = density_scaled_k(n_localities, fraction=0.15)
    with geometry_csv.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDS)
        writer.writeheader()
        for species_index, species in enumerate(species_names):
            for locality in range(n_localities):
                writer.writerow(
                    {
                        "species": species,
                        "locality_index": locality,
                        "x_km": species_index * 1000.0 + locality * 11.0,
                        "y_km": float((locality * 7) % 13),
                        "z_km": float((locality * 5) % 17),
                        "graph_k": graph_k,
                    }
                )

    table = load_frozen_genetic_geometry_csv(geometry_csv)
    fingerprint = geometry_fingerprint(
        [
            SpeciesGeometry(species=name, coordinates=table.geometries[name].coordinates)
            for name in table.species
        ]
    )
    train = species_names[:75]
    evaluation = species_names[75:]
    phase2 = {
        "schema": "ttf_genetic_phylogatr_confirmatory_phase2_mask_v0.1",
        "status": "PASS_TO_SYNTHETIC_GATE",
        "phase1": {
            "dataset_digest_sha256": "a" * 64,
        },
        "mask_contract": {
            "nucleotide_identity_persisted": False,
            "pairwise_nucleotide_differences_computed": False,
        },
        "species": {
            "survivors": 150,
        },
        "split": {
            "inherit_phase1_without_resplitting": True,
            "train_species": list(train),
            "eval_species": list(evaluation),
        },
        "geometry_csv_sha256": sha256_path(geometry_csv),
        "geometry_fingerprint_sha256": fingerprint,
        "character_mask_opened": True,
        "confirmatory_sequence_identity_opened": False,
        "confirmatory_pairwise_genetic_distances_opened": False,
        "confirmatory_ttf_statistic_opened": False,
    }
    phase2_path = tmp_path / "phase2.json"
    phase2_path.write_text(json.dumps(phase2, indent=2, sort_keys=True) + "\n")
    return geometry_csv, phase2_path


def test_retention_rule_is_deterministic_nested_and_respects_floor() -> None:
    digest = "d" * 64
    fractions = [1.0, 0.875, 0.75, 0.625, 0.5]
    selections = [
        retained_locality_indices(
            digest,
            "Freshus alpha",
            24,
            fraction,
            minimum_localities=12,
        )
        for fraction in fractions
    ]
    assert [len(selection) for selection in selections] == [24, 21, 18, 15, 12]
    assert selections == [
        retained_locality_indices(
            digest,
            "Freshus alpha",
            24,
            fraction,
            minimum_localities=12,
        )
        for fraction in fractions
    ]
    for higher, lower in zip(selections, selections[1:]):
        assert set(lower).issubset(higher)


def test_fragility_plan_preserves_formal_panel_and_has_no_authority(tmp_path: Path) -> None:
    geometry_csv, phase2_path = _fixture(tmp_path)
    plan = build_fragility_plan(
        geometry_csv,
        phase2_path,
        PHASE3_RULE,
        FRAGILITY_RULE,
    )
    phase2 = json.loads(phase2_path.read_text())

    assert [level.retention_fraction for level in plan.levels] == [
        1.0,
        0.875,
        0.75,
        0.625,
        0.5,
    ]
    assert plan.levels[0].geometry_fingerprint_sha256 == phase2[
        "geometry_fingerprint_sha256"
    ]
    assert all(level.metrics["species_count"] == 150 for level in plan.levels)
    assert all(level.metrics["minimum_localities"] >= 12 for level in plan.levels)
    assert set(plan.train_species) == set(phase2["split"]["train_species"])
    assert set(plan.eval_species) == set(phase2["split"]["eval_species"])

    receipt_path = write_fragility_plan(plan, tmp_path / "fragility")
    receipt = json.loads(receipt_path.read_text())
    assert receipt["formal_gate_d_decision_made_by_this_receipt"] is False
    assert receipt["phase4_identity_opening_authorized_by_this_receipt"] is False
    assert all(value is False for value in receipt["authority_firewall"].values())
    assert receipt["confirmatory_sequence_identity_opened"] is False
    assert receipt["confirmatory_pairwise_genetic_distances_opened"] is False
    assert receipt["confirmatory_ttf_statistic_opened"] is False

    for higher, lower in zip(receipt["levels"], receipt["levels"][1:]):
        for species in receipt["species"]["train_species"] + receipt["species"]["eval_species"]:
            assert set(lower["source_locality_indices"][species]).issubset(
                higher["source_locality_indices"][species]
            )


def test_fragility_plan_rejects_any_fresh_identity_opening(tmp_path: Path) -> None:
    geometry_csv, phase2_path = _fixture(tmp_path)
    phase2 = json.loads(phase2_path.read_text())
    phase2["confirmatory_sequence_identity_opened"] = True
    phase2_path.write_text(json.dumps(phase2, indent=2, sort_keys=True) + "\n")

    with pytest.raises(RuntimeError, match="forbidden after fresh nucleotide identity opening"):
        build_fragility_plan(
            geometry_csv,
            phase2_path,
            PHASE3_RULE,
            FRAGILITY_RULE,
        )
