#!/usr/bin/env python3
from __future__ import annotations

import argparse
from dataclasses import asdict
import json
from pathlib import Path

import numpy as np

from ttf.calibration import CalibrationCell, seed_for
from ttf.controls import paired_opportunity_control_test
from ttf.core import split_species
from ttf.nulls import edges_on_fixed_graphs, fixed_graphs
from ttf.nuisance import residualize_edge_sets
from ttf.semisynthetic import simulate_on_geometry

EARTH_RADIUS_KM = 6371.0088
VARIANTS = {"paired_raw", "paired_ibd"}


def load_geometry(path: Path) -> dict[str, np.ndarray]:
    payload = json.loads(path.read_text())
    if payload.get("schema") != "ttf_iran_plateau_geometry_v0.1":
        raise RuntimeError("unexpected Iranian geometry schema")
    raw = {
        str(item["species"]): np.asarray(
            [[row["latitude"], row["longitude"]] for row in item["populations"]],
            dtype=float,
        )
        for item in payload["species"]
    }
    pooled = np.vstack(tuple(raw.values()))
    mean_lat = float(pooled[:, 0].mean())
    mean_lon = float(pooled[:, 1].mean())
    lat0 = float(np.deg2rad(mean_lat))
    out = {}
    for species, latlon in raw.items():
        lat = np.deg2rad(latlon[:, 0])
        lon = np.deg2rad(latlon[:, 1])
        x = EARTH_RADIUS_KM * np.cos(lat0) * (lon - np.deg2rad(mean_lon))
        y = EARTH_RADIUS_KM * (lat - np.deg2rad(mean_lat))
        out[species] = np.column_stack([x, y])
    return out


def main() -> int:
    parser = argparse.ArgumentParser(description="Calibrate opportunity-incremental TTF on Iranian geometry.")
    parser.add_argument("--geometry", type=Path, required=True)
    parser.add_argument("--variant", choices=sorted(VARIANTS), required=True)
    parser.add_argument("--shared-fraction", type=float, required=True)
    parser.add_argument("--amplitude", type=float, required=True)
    parser.add_argument("--replicates", type=int, default=100)
    parser.add_argument("--bootstrap", type=int, default=1999)
    parser.add_argument("--alpha", type=float, default=0.05)
    parser.add_argument("--k", type=int, default=4)
    parser.add_argument("--noise-sd", type=float, default=0.8)
    parser.add_argument("--transition-width", type=float, default=0.20)
    parser.add_argument("--min-crossing-species", type=int, default=6)
    parser.add_argument("--split-seed", type=int, default=20260907)
    parser.add_argument("--ibd-degree", type=int, default=1)
    parser.add_argument("--seed", type=int, default=20260907)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    geometry = load_geometry(args.geometry)
    labels = tuple(sorted(geometry))
    train, evaluation = split_species(labels, eval_fraction=6.0 / len(labels), seed=args.split_seed)
    if len(train) != 3 or len(evaluation) != 6:
        raise RuntimeError("paired-control benchmark requires the frozen 3/6 split")

    p_values: list[float] = []
    deltas: list[float] = []
    boundary_stats: list[float] = []
    opportunity_stats: list[float] = []
    bandwidths: list[float] = []
    for replicate in range(args.replicates):
        # Reuse the exact world seed family from the prior ablation.  Only the
        # statistic changes; the synthetic worlds do not.
        world = simulate_on_geometry(
            geometry,
            shared_fraction=args.shared_fraction,
            amplitude=args.amplitude,
            noise_sd=args.noise_sd,
            transition_width=args.transition_width,
            min_shared_crossing_species=args.min_crossing_species,
            k=args.k,
            seed=seed_for(args.seed, "iran-ablation-world", args.shared_fraction, args.amplitude, replicate),
        )
        graphs = fixed_graphs(world.samples, k=args.k)
        edge_map = edges_on_fixed_graphs(world.samples, graphs)
        edge_sets = tuple(edge_map[name] for name in labels)
        if args.variant == "paired_ibd":
            edge_sets = residualize_edge_sets(edge_sets, degree=args.ibd_degree)
        edge_map = {edges.species: edges for edges in edge_sets}

        result = paired_opportunity_control_test(
            [edge_map[name] for name in train],
            [edge_map[name] for name in evaluation],
            bandwidth=world.bandwidth,
            n_bootstrap=args.bootstrap,
            seed=seed_for(
                args.seed,
                "iran-paired-bootstrap",
                args.variant,
                args.shared_fraction,
                args.amplitude,
                replicate,
            ),
        )
        p_values.append(result.p_value)
        deltas.append(result.observed.statistic)
        boundary_stats.append(float(np.mean([
            value for value in result.observed.boundary_scores.values() if np.isfinite(value)
        ])))
        opportunity_stats.append(float(np.mean([
            value for value in result.observed.opportunity_scores.values() if np.isfinite(value)
        ])))
        bandwidths.append(world.bandwidth)

    if not np.allclose(bandwidths, bandwidths[0], atol=1e-12, rtol=0.0):
        raise RuntimeError("geometry bandwidth drifted")
    p = np.asarray(p_values, dtype=float)
    cell = CalibrationCell(
        shared_fraction=float(args.shared_fraction),
        amplitude=float(args.amplitude),
        n_replicates=int(args.replicates),
        alpha=float(args.alpha),
        rejection_rate=float(np.mean(p <= args.alpha)),
        mean_statistic=float(np.mean(deltas)),
        mean_null_statistic=0.0,
        median_p_value=float(np.median(p)),
    )
    payload = {
        "schema": "ttf_iran_paired_control_cell_v0.1",
        "variant": args.variant,
        "cell": asdict(cell),
        "diagnostic_means": {
            "boundary_transfer": float(np.mean(boundary_stats)),
            "opportunity_transfer": float(np.mean(opportunity_stats)),
            "increment": float(np.mean(deltas)),
        },
        "design": {
            "statistic": "mean_s[Spearman(boundary_exposure, turnover) - Spearman(opportunity_exposure, turnover)]",
            "opportunity_control_trait_free": True,
            "geometry_schema": "ttf_iran_plateau_geometry_v0.1",
            "bandwidth": float(bandwidths[0]),
            "bandwidth_rule": "pooled median species-local kNN edge length after normalization",
            "k": int(args.k),
            "noise_sd": float(args.noise_sd),
            "transition_width": float(args.transition_width),
            "min_shared_crossing_species": int(args.min_crossing_species),
            "bootstrap_resamples": int(args.bootstrap),
            "train_species": list(train),
            "eval_species": list(evaluation),
            "ibd_residualization": "leave-one-edge-out rank-linear edge-length nuisance" if args.variant == "paired_ibd" else None,
            "ibd_degree": int(args.ibd_degree) if args.variant == "paired_ibd" else None,
            "world_seed_family": "iran-ablation-world",
            "worlds_paired_to_previous_ablation": True,
            "master_seed": int(args.seed),
            "empirical_genetic_outcomes_used": False,
            "known_boundary_labels_used": False,
        },
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"variant": args.variant, "cell": payload["cell"], "diagnostic_means": payload["diagnostic_means"]}, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
