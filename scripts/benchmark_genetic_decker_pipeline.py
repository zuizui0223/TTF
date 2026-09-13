#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import time
from collections import defaultdict
from pathlib import Path

import numpy as np

from ttf.genetic_gate import prepare_genetic_ttf_design, score_genetic_world
from ttf.genetic_geometry import GeneticSamplingGeometry, prepare_density_scaled_genetic_geometry
from ttf.genetic_simulate import simulate_genetic_distance_world
from ttf.geometry import SpeciesGeometry, geometry_fingerprint


def sha256_path(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def load_geometries(path: Path) -> dict[str, GeneticSamplingGeometry]:
    rows: dict[str, list[tuple[int, np.ndarray, int]]] = defaultdict(list)
    with path.open(newline="") as handle:
        for row in csv.DictReader(handle):
            species = str(row["species"])
            rows[species].append(
                (
                    int(row["locality_index"]),
                    np.asarray(
                        [float(row["x_km"]), float(row["y_km"]), float(row["z_km"])],
                        dtype=float,
                    ),
                    int(row["graph_k"]),
                )
            )
    out: dict[str, GeneticSamplingGeometry] = {}
    for species in sorted(rows):
        ordered = sorted(rows[species], key=lambda item: item[0])
        if [item[0] for item in ordered] != list(range(len(ordered))):
            raise RuntimeError(f"locality index drift for {species}")
        coordinates = np.vstack([item[1] for item in ordered])
        k_values = {item[2] for item in ordered}
        if len(k_values) != 1:
            raise RuntimeError(f"graph k drift within {species}")
        geometry = prepare_density_scaled_genetic_geometry(coordinates, neighbor_fraction=0.15)
        if geometry.graph_k != next(iter(k_values)):
            raise RuntimeError(f"density-scaled k reconstruction drift for {species}")
        out[species] = geometry
    return out


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
    if authorization.get("schema") != "ttf_genetic_decker_execution_benchmark_authorization_v0.1":
        raise RuntimeError("benchmark authorization drift")
    if authorization.get("status") != "authorize_exactly_one_synthetic_execution_benchmark":
        raise RuntimeError("benchmark is not authorized")
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

    cases = (
        ("ibd_only", 0.0, 0.0, 0.0, 2026091301),
        ("private_A2", 0.0, 2.0, 0.10, 2026091302),
        ("shared_A2", 1.0, 2.0, 0.10, 2026091303),
    )
    case_results = []
    for label, shared, amplitude, noise_sd, seed in cases:
        world = simulate_genetic_distance_world(
            geometries,
            shared_fraction=shared,
            residual_amplitude=amplitude,
            ibd_strength=1.0,
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

    ibd_only = next(row for row in case_results if row["label"] == "ibd_only")
    if abs(float(ibd_only["statistic"])) > 1e-12:
        raise RuntimeError("IBD-only world did not collapse to zero held-out transfer")
    if abs(float(ibd_only["training_strength"])) > 1e-12:
        raise RuntimeError("IBD-only world did not collapse to zero training strength")

    output = {
        "schema": "ttf_genetic_decker_execution_benchmark_v0.1",
        "status": "synthetic_execution_benchmark_only_not_qualification",
        "authorization": str(args.authorization),
        "geometry_fingerprint_sha256": fingerprint,
        "species_count": len(geometries),
        "train_species": len(design.train_species),
        "eval_species": len(design.eval_species),
        "total_edges": int(sum(item.n_edges for item in design.template_edges.values())),
        "design_preparation_seconds": float(design_seconds),
        "cases": case_results,
        "empirical_genetic_outcomes_opened": False,
        "qualification_claim_made": False,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(output, indent=2, sort_keys=True) + "\n")
    print(json.dumps(output, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
