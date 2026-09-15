from __future__ import annotations

import csv
import json
from pathlib import Path

import pytest

from ttf.genetic_geometry_io import load_frozen_genetic_geometry_csv, sha256_path
from ttf.geometry import SpeciesGeometry, geometry_fingerprint
from ttf.heterogeneous_inference import density_scaled_k
from ttf.phylogatr_fragility import build_fragility_plan, write_fragility_plan
from ttf.phylogatr_fragility_execution import (
    build_fragility_execution_plan,
    derive_fragility_level_master_seed,
    load_fragility_execution_context,
)


PHASE3_RULE = Path("docs/supporting/genetic_phylogatr_phase3_gate_d_rule_v0.1.json")
FRAGILITY_RULE = Path(
    "docs/supporting/genetic_phylogatr_gate_d_fragility_diagnostic_v0.1.json"
)
EXECUTION_RULE = Path(
    "docs/supporting/genetic_phylogatr_gate_d_fragility_execution_v0.1.json"
)
FIELDS = ["species", "locality_index", "x_km", "y_km", "z_km", "graph_k"]


def _phase2_fixture(tmp_path: Path) -> tuple[Path, Path]:
    geometry_csv = tmp_path / "phase2_geometry.csv"
    species_names = tuple(f"fresh_exec_{index:03d}" for index in range(150))
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
                        "x_km": species_index * 1000.0 + locality * 13.0,
                        "y_km": float((locality * 11) % 19),
                        "z_km": float((locality * 7) % 23),
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
        "phase1": {"dataset_digest_sha256": "a" * 64},
        "mask_contract": {
            "nucleotide_identity_persisted": False,
            "pairwise_nucleotide_differences_computed": False,
        },
        "species": {"survivors": 150},
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


def _geometry_plan(tmp_path: Path) -> tuple[Path, Path, Path]:
    geometry_csv, phase2_path = _phase2_fixture(tmp_path)
    plan = build_fragility_plan(
        geometry_csv,
        phase2_path,
        PHASE3_RULE,
        FRAGILITY_RULE,
    )
    output_dir = tmp_path / "geometry_plan"
    receipt = write_fragility_plan(plan, output_dir)
    return output_dir, receipt, phase2_path


def _formal_qualification(tmp_path: Path, phase2_path: Path) -> Path:
    phase2 = json.loads(phase2_path.read_text())
    cells = []
    for shared, amplitude in [(0.0, 0.5), (0.0, 1.0), (0.0, 2.0), (0.0, 3.0), (1.0, 2.0)]:
        if shared == 0.0 and amplitude == 3.0:
            rate, low, high = 0.06, 0.04, 0.10033475332223055
        elif shared == 1.0 and amplitude == 2.0:
            rate, low, high = 0.84, 0.805, 0.87
        else:
            rate, low, high = 0.05, 0.03, 0.08
        cells.append(
            {
                "shared_fraction": shared,
                "residual_amplitude": amplitude,
                "n_worlds": 500,
                "rejections": int(round(rate * 500)),
                "rejection_rate": rate,
                "wilson95_low": low,
                "wilson95_high": high,
                "selected_pair_counts": {},
            }
        )
    payload = {
        "schema": "ttf_genetic_phylogatr_phase3_qualification_v0.1",
        "status": "NOT_EVALUABLE",
        "geometry_fingerprint_sha256": phase2["geometry_fingerprint_sha256"],
        "cells": cells,
        "type1_gate": {
            "wilson95_upper_ceiling": 0.10,
            "max_observed_wilson95_upper": 0.10033475332223055,
            "pass": False,
        },
        "power_gate": {
            "cell": {"shared_fraction": 1.0, "residual_amplitude": 2.0},
            "wilson95_lower_floor": 0.80,
            "observed_wilson95_lower": 0.805,
            "pass": True,
        },
        "passed": False,
        "failure_interpretation": "NOT_EVALUABLE_under_this_fresh_geometry_not_absence_of_transferable_phylogeography",
        "phase4_identity_opening_eligible": False,
        "confirmatory_sequence_identity_opened": False,
        "confirmatory_pairwise_genetic_distances_opened": False,
        "confirmatory_ttf_statistic_opened": False,
    }
    path = tmp_path / "formal_qualification.json"
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    return path


def test_execution_plan_imports_formal_anchor_and_never_reruns_retention_one(
    tmp_path: Path,
) -> None:
    geometry_dir, receipt, phase2_path = _geometry_plan(tmp_path)
    formal = _formal_qualification(tmp_path, phase2_path)
    plan = build_fragility_execution_plan(
        receipt,
        phase2_path,
        PHASE3_RULE,
        formal,
        EXECUTION_RULE,
    )
    assert plan["formal_anchor"]["retention_fraction"] == 1.0
    assert plan["formal_anchor"]["formal_status"] == "NOT_EVALUABLE"
    assert plan["formal_anchor"]["A3_private_null_Wilson95_upper"] == pytest.approx(
        0.10033475332223055
    )
    assert [level["retention_fraction"] for level in plan["levels"]] == [
        0.875,
        0.75,
        0.625,
        0.5,
    ]
    assert all(level["retention_fraction"] != 1.0 for level in plan["levels"])
    assert all(level["formal_gate_d_decision_authority"] is False for level in plan["levels"])
    assert plan["formal_gate_d_decision_made_by_this_plan"] is False
    assert plan["phase4_identity_opening_authorized_by_this_plan"] is False
    assert all(value is False for value in plan["authority_firewall"].values())

    runnable = [level for level in plan["levels"] if level["status"] == "SYNTHETIC_DIAGNOSTIC_RUNNABLE"]
    assert runnable
    assert all(level["reference_job_count"] == 64 for level in runnable)
    assert all(level["observed_job_count"] == 10 for level in runnable)
    assert len({level["master_seed"] for level in plan["levels"]}) == 4

    plan_path = tmp_path / "execution_plan.json"
    plan_path.write_text(json.dumps(plan, indent=2, sort_keys=True) + "\n")
    context = load_fragility_execution_context(
        plan_path,
        geometry_dir,
        phase2_path,
        PHASE3_RULE,
        runnable[0]["retention_fraction"],
    )
    assert context.retention_fraction == runnable[0]["retention_fraction"]
    assert len(context.train_species) == 75
    assert len(context.eval_species) == 75


def test_fragility_level_seed_is_deterministic_and_retention_specific() -> None:
    seed_a = derive_fragility_level_master_seed("a" * 64, "b" * 64, 0.875, "c" * 64)
    seed_b = derive_fragility_level_master_seed("a" * 64, "b" * 64, 0.875, "c" * 64)
    seed_c = derive_fragility_level_master_seed("a" * 64, "b" * 64, 0.75, "c" * 64)
    assert seed_a == seed_b
    assert seed_a != seed_c


def test_execution_plan_rejects_formal_receipt_after_identity_opening(tmp_path: Path) -> None:
    _, receipt, phase2_path = _geometry_plan(tmp_path)
    formal = _formal_qualification(tmp_path, phase2_path)
    payload = json.loads(formal.read_text())
    payload["confirmatory_sequence_identity_opened"] = True
    formal.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    with pytest.raises(RuntimeError, match="formal qualification firewall is open"):
        build_fragility_execution_plan(
            receipt,
            phase2_path,
            PHASE3_RULE,
            formal,
            EXECUTION_RULE,
        )
