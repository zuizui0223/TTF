#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

from ttf.genetic_geometry_io import load_frozen_genetic_geometry_csv, sha256_path
from ttf.geometry import SpeciesGeometry, geometry_fingerprint
from ttf.phylogatr_phase3 import derive_phase3_master_seed


CRITICAL_CODE_PATHS = (
    "src/ttf/genetic_simulate.py",
    "src/ttf/genetic_ibd.py",
    "src/ttf/genetic_gate.py",
    "src/ttf/genetic_geometry.py",
    "src/ttf/genetic_geometry_io.py",
    "src/ttf/genetic_batch_execution.py",
    "src/ttf/cached_chunked_transfer.py",
    "src/ttf/chunked_transfer.py",
    "src/ttf/batch.py",
    "src/ttf/geometry_control.py",
    "src/ttf/private_strength.py",
    "src/ttf/profiled_private_null.py",
    "src/ttf/precision.py",
    "src/ttf/calibration.py",
    "src/ttf/phylogatr_phase3.py",
    "scripts/authorize_phylogatr_phase3_gate_d.py",
    "scripts/plan_phylogatr_phase3_shards.py",
    "scripts/run_phylogatr_phase3_reference_shard.py",
    "scripts/aggregate_phylogatr_phase3_references.py",
    "scripts/run_phylogatr_phase3_observed_shard.py",
    "scripts/aggregate_phylogatr_phase3_qualification.py",
)


def _load_json(path: Path, schema: str) -> dict:
    payload = json.loads(path.read_text())
    if payload.get("schema") != schema:
        raise RuntimeError(f"unexpected schema for {path}: {payload.get('schema')!r}")
    return payload


def _assert_false_mapping(payload: dict, label: str) -> None:
    if not payload or any(value is not False for value in payload.values()):
        raise RuntimeError(f"{label} is not fully closed")


def _assert_phase2_blindness(phase2: dict) -> None:
    if phase2.get("character_mask_opened") is not True:
        raise RuntimeError("Phase-2 character mask was not opened/frozen")
    if phase2.get("confirmatory_sequence_identity_opened") is not False:
        raise RuntimeError("cannot authorize Phase 3 after fresh nucleotide identity opening")
    if phase2.get("confirmatory_pairwise_genetic_distances_opened") is not False:
        raise RuntimeError("cannot authorize Phase 3 after fresh genetic-distance opening")
    if phase2.get("confirmatory_ttf_statistic_opened") is not False:
        raise RuntimeError("cannot authorize Phase 3 after fresh empirical TTF opening")
    mask = phase2.get("mask_contract")
    if not isinstance(mask, dict):
        raise RuntimeError("Phase-2 mask contract missing")
    if mask.get("nucleotide_identity_persisted") is not False:
        raise RuntimeError("Phase-2 mask contract persisted nucleotide identity")
    if mask.get("pairwise_nucleotide_differences_computed") is not False:
        raise RuntimeError("Phase-2 mask contract computed nucleotide differences")
    split = phase2.get("split")
    if not isinstance(split, dict) or split.get("inherit_phase1_without_resplitting") is not True:
        raise RuntimeError("Phase-2 split was not inherited unchanged from Phase 1")


def main() -> int:
    ap = argparse.ArgumentParser(
        description="Authorize fresh phylogatR Phase-3 Gate-D on one exact Phase-2 survivor geometry."
    )
    ap.add_argument("--geometry", type=Path, required=True)
    ap.add_argument("--phase2-manifest", type=Path, required=True)
    ap.add_argument("--phase3-rule", type=Path, required=True)
    ap.add_argument("--repo-root", type=Path, default=Path("."))
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()

    phase2 = _load_json(
        args.phase2_manifest,
        "ttf_genetic_phylogatr_confirmatory_phase2_mask_v0.1",
    )
    rule = _load_json(
        args.phase3_rule,
        "ttf_genetic_phylogatr_phase3_gate_d_rule_v0.1",
    )
    required_status = str(rule["geometry_contract"]["required_phase2_status"])
    if phase2.get("status") != required_status:
        raise RuntimeError(
            f"fresh Phase-3 authorization requires {required_status}, got {phase2.get('status')!r}"
        )
    _assert_phase2_blindness(phase2)
    _assert_false_mapping(rule["outcome_firewall"], "phase-3 rule outcome firewall")

    geometry_sha = sha256_path(args.geometry)
    if geometry_sha != phase2["geometry_csv_sha256"]:
        raise RuntimeError("Phase-2 geometry CSV digest drift before authorization")
    table = load_frozen_genetic_geometry_csv(
        args.geometry,
        expected_sha256=geometry_sha,
        expected_neighbor_fraction=float(rule["geometry_contract"]["neighbor_fraction"]),
    )
    minimum_training_edges = int(
        rule["geometry_contract"]["minimum_endpoint_disjoint_ibd_training_edges"]
    )
    insufficient = [
        name
        for name, geometry in table.geometries.items()
        if geometry.min_endpoint_disjoint_training_edges < minimum_training_edges
    ]
    if insufficient:
        raise RuntimeError(
            "fresh Phase-2 geometry lost endpoint-disjoint IBD support: "
            + ", ".join(sorted(insufficient)[:10])
        )

    train = tuple(map(str, phase2["split"]["train_species"]))
    evaluation = tuple(map(str, phase2["split"]["eval_species"]))
    species = set(table.species)
    if not train or not evaluation or set(train) & set(evaluation):
        raise RuntimeError("fresh Phase-2 split is empty or overlapping")
    if set(train) | set(evaluation) != species:
        raise RuntimeError("fresh Phase-2 split does not cover the survivor geometry exactly")
    if len(species) != int(phase2["species"]["survivors"]):
        raise RuntimeError("fresh Phase-2 survivor count drift")
    minimum = int(rule["geometry_contract"]["minimum_surviving_species"])
    if len(species) < minimum:
        raise RuntimeError("fresh Phase-2 panel is below the frozen Phase-3 minimum")

    fingerprint = geometry_fingerprint(
        [
            SpeciesGeometry(species=name, coordinates=table.geometries[name].coordinates)
            for name in sorted(table.geometries)
        ]
    )
    if fingerprint != phase2["geometry_fingerprint_sha256"]:
        raise RuntimeError("fresh Phase-2 geometry fingerprint drift before authorization")

    dataset_digest = str(phase2["phase1"]["dataset_digest_sha256"])
    master_seed = derive_phase3_master_seed(dataset_digest, fingerprint)
    execution = dict(rule["execution_defaults"])
    qualification = rule["qualification"]
    configs = rule["synthetic_world"]["private_reference_configurations"]
    reference_n = int(qualification["reference_worlds_per_configuration"])
    observed_n = int(qualification["observed_worlds_per_cell"])
    reference_shard = int(execution["reference_shard_size"])
    observed_shard = int(execution["observed_shard_size"])
    reference_jobs = len(configs) * math.ceil(reference_n / reference_shard)
    observed_jobs = len(qualification["mandatory_primary_cells"]) * math.ceil(
        observed_n / observed_shard
    )

    repo_root = args.repo_root.resolve()
    frozen_code: dict[str, str] = {}
    for relative in CRITICAL_CODE_PATHS:
        path = repo_root / relative
        if not path.is_file():
            raise RuntimeError(f"critical fresh Phase-3 code file missing: {relative}")
        frozen_code[relative] = sha256_path(path)

    out = {
        "schema": "ttf_genetic_phylogatr_phase3_gate_d_authorization_v0.1",
        "status": "authorize_frozen_fresh_phylogatr_phase3_gate_d",
        "purpose": "Authorize synthetic Gate-D only on this exact response-blind Phase-2 survivor geometry before fresh nucleotide identity is opened.",
        "phase2_manifest_sha256": sha256_path(args.phase2_manifest),
        "phase3_rule_sha256": sha256_path(args.phase3_rule),
        "phase2_mask_rule_sha256": phase2["rule_sha256"],
        "dataset_digest_sha256": dataset_digest,
        "geometry_csv_sha256": geometry_sha,
        "geometry_fingerprint_sha256": fingerprint,
        "species": {
            "survivors": len(species),
            "train": len(train),
            "eval": len(evaluation),
            "train_species": list(train),
            "eval_species": list(evaluation),
        },
        "master_seed": master_seed,
        "seed_rule": rule["seed_contract"],
        "execution": {
            "batch_width": int(execution["batch_width"]),
            "edge_chunk_size": int(execution["edge_chunk_size"]),
            "train_chunk_size": int(execution["train_chunk_size"]),
            "reference_shard_size": reference_shard,
            "observed_shard_size": observed_shard,
            "reference_expected_jobs": reference_jobs,
            "observed_expected_jobs": observed_jobs,
            "sharding_semantics": execution["sharding_semantics"],
        },
        "frozen_inference": {
            "private_reference_configurations": list(configs.keys()),
            "reference_worlds_per_configuration": reference_n,
            "profile_strength_draws": int(qualification["profile_strength_draws"]),
            "calibration_statistic_draws": int(qualification["calibration_statistic_draws"]),
            "selected_private_configurations": int(
                rule["core_method"]["profiled_private_selected_configurations"]
            ),
            "mandatory_primary_cells": qualification["mandatory_primary_cells"],
            "observed_worlds_per_cell": observed_n,
            "alpha": float(qualification["alpha"]),
            "type1_wilson95_upper_ceiling": float(
                qualification["type1_wilson95_upper_ceiling"]
            ),
            "shared_A2_wilson95_lower_floor": float(
                qualification["shared_A2_wilson95_lower_floor"]
            ),
            "failure_interpretation": qualification["failure_interpretation"],
        },
        "frozen_code_sha256": frozen_code,
        "outcome_firewall": dict(rule["outcome_firewall"]),
        "forbidden_after_authorization": [
            "changing the Phase-2 survivor species set, graph, or inherited split",
            "changing the Phase-3 rule, master-seed derivation, world counts, reference family, alpha, or Wilson gates",
            "changing the bandwidth or IBD residualization based on synthetic qualification results",
            "opening fresh nucleotide identity, pairwise genetic distances, or empirical TTF values before a complete PASS receipt",
        ],
        "claim_boundary": "This authorization permits synthetic qualification only. It does not open nucleotide identity and does not constitute an empirical genetic result.",
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(out, indent=2, sort_keys=True) + "\n")
    print(
        json.dumps(
            {
                "status": out["status"],
                "geometry_fingerprint_sha256": fingerprint,
                "species": len(species),
                "train": len(train),
                "eval": len(evaluation),
                "master_seed": master_seed,
                "reference_expected_jobs": reference_jobs,
                "observed_expected_jobs": observed_jobs,
                "fresh_nucleotide_identity_opened": False,
            },
            sort_keys=True,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
