#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import time
from pathlib import Path

import numpy as np

from benchmark_genetic_decker_pipeline import load_geometries, sha256_path
from ttf.genetic_gate import prepare_genetic_ttf_design, score_genetic_world
from ttf.genetic_simulate import simulate_genetic_distance_world
from ttf.geometry import SpeciesGeometry, geometry_fingerprint


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--geometry", type=Path, required=True)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--source-ledger", type=Path, required=True)
    parser.add_argument("--rule", type=Path, required=True)
    parser.add_argument("--authorization", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    manifest = json.loads(args.manifest.read_text())
    ledger = json.loads(args.source_ledger.read_text())
    rule = json.loads(args.rule.read_text())
    authorization = json.loads(args.authorization.read_text())

    if manifest.get("schema") != "ttf_genetic_decker_pilot_geometry_v0.2":
        raise RuntimeError("canonical v0.2 manifest required")
    if manifest.get("synthetic_worlds_scored_before_freeze") != 0:
        raise RuntimeError("geometry was performance-opened before freeze")
    if ledger.get("schema") != "ttf_genetic_decker_pilot_geometry_source_v0.2":
        raise RuntimeError("canonical v0.2 source ledger required")
    if authorization.get("schema") != "ttf_genetic_decker_execution_benchmark_authorization_v0.2":
        raise RuntimeError("benchmark authorization drift")
    if authorization.get("status") != "authorize_exactly_one_v02_synthetic_execution_benchmark":
        raise RuntimeError("benchmark is not authorized")
    if authorization.get("geometry_fingerprint_sha256") != ledger.get("geometry_fingerprint_sha256"):
        raise RuntimeError("authorization/ledger fingerprint drift")
    if authorization.get("geometry_csv_sha256") != ledger.get("geometry_csv_sha256"):
        raise RuntimeError("authorization/ledger CSV drift")
    if ledger.get("geometry_csv_sha256") != sha256_path(args.geometry):
        raise RuntimeError("geometry CSV differs from frozen source ledger")
    if ledger.get("geometry_fingerprint_sha256") != manifest.get("geometry_fingerprint_sha256"):
        raise RuntimeError("manifest/ledger fingerprint drift")
    if ledger.get("rule_sha256") != sha256_path(args.rule):
        raise RuntimeError("rule differs from frozen source ledger")

    geometries = load_geometries(args.geometry)
    fingerprint = geometry_fingerprint(
        [SpeciesGeometry(species=name, coordinates=geometries[name].coordinates) for name in sorted(geometries)]
    )
    if fingerprint != ledger["geometry_fingerprint_sha256"]:
        raise RuntimeError("reconstructed geometry fingerprint drift")

    split = manifest["split"]
    design_started = time.perf_counter()
    design = prepare_genetic_ttf_design(
        geometries,
        train_species=split["train_species"],
        eval_species=split["eval_species"],
        bandwidth=float(rule["core_v011_inheritance"]["bandwidth_km"]),
        prior_strength=float(rule["core_v011_inheritance"]["prior_strength"]),
        segment_points=int(rule["core_v011_inheritance"]["segment_points"]),
        min_training_edges=int(rule["geometry_filter"]["min_endpoint_disjoint_ibd_training_edges"]),
        strength_neighbours=int(rule["core_v011_inheritance"]["strength_neighbours"]),
    )
    design_seconds = time.perf_counter() - design_started

    cases = tuple(
        (
            str(case["label"]),
            float(case["shared_fraction"]),
            float(case["residual_amplitude"]),
            float(case["noise_sd"]),
            int(case["seed"]),
        )
        for case in authorization["authorized_cases"]
    )
    case_results = []
    for label, shared, amplitude, noise_sd, seed in cases:
        world = simulate_genetic_distance_world(
            geometries,
            shared_fraction=shared,
            residual_amplitude=amplitude,
            ibd_strength=float(rule["genetic_world"]["ibd_strength_primary"]),
            noise_sd=noise_sd,
            transition_width=float(rule["genetic_world"]["transition_width"]),
            noise_dimensions=int(rule["genetic_world"]["noise_dimensions"]),
            seed=seed,
        )
        started = time.perf_counter()
        score = score_genetic_world(design, world)
        elapsed = time.perf_counter() - started
        finite_scores = [value for value in score.species_scores.values() if np.isfinite(value)]
        case_results.append(
            {
                "label": label,
                "shared_fraction": shared,
                "residual_amplitude": amplitude,
                "noise_sd": noise_sd,
                "statistic": float(score.statistic),
                "training_strength": float(score.training_strength),
                "finite_eval_species_scores": len(finite_scores),
                "score_seconds": float(elapsed),
            }
        )

    invariants = authorization["required_invariants"]
    ibd_only = next(row for row in case_results if row["label"] == "ibd_only")
    if abs(float(ibd_only["statistic"])) > float(invariants["ibd_only_statistic_absolute_max"]):
        raise RuntimeError("IBD-only world did not collapse to zero held-out transfer")
    if abs(float(ibd_only["training_strength"])) > float(invariants["ibd_only_training_strength_absolute_max"]):
        raise RuntimeError("IBD-only world did not collapse to zero training strength")
    if invariants["all_110_eval_species_scores_finite"] and any(
        int(row["finite_eval_species_scores"]) != 110 for row in case_results
    ):
        raise RuntimeError("not all held-out species received finite scores")

    output = {
        "schema": "ttf_genetic_decker_execution_benchmark_v0.2",
        "status": "synthetic_execution_benchmark_only_not_qualification",
        "authorization": str(args.authorization),
        "source_ledger": str(args.source_ledger),
        "geometry_fingerprint_sha256": fingerprint,
        "species_count": len(geometries),
        "train_species": len(design.train_species),
        "eval_species": len(design.eval_species),
        "total_edges": int(sum(item.n_edges for item in design.template_edges.values())),
        "design_preparation_seconds": float(design_seconds),
        "cases": case_results,
        "empirical_genetic_outcomes_opened": False,
        "qualification_claim_made": False,
        "v0.1_execution_benchmark_remains_invalid": True,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(output, indent=2, sort_keys=True) + "\n")
    print(json.dumps(output, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
