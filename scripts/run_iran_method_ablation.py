#!/usr/bin/env python3
from __future__ import annotations

import argparse
from dataclasses import asdict
import json
from pathlib import Path

import numpy as np

from ttf.calibration import CalibrationCell, seed_for
from ttf.core import split_species
from ttf.crossfit import balanced_species_folds, crossfit_species_bootstrap_test
from ttf.inference import SpeciesBootstrapResult, centered_species_bootstrap_mean_test
from ttf.nulls import edges_on_fixed_graphs, fixed_graphs
from ttf.nuisance import residualize_edge_sets
from ttf.semisynthetic import simulate_on_geometry
from ttf.transfer import TransferResult, prepare_transfer

EARTH_RADIUS_KM = 6371.0088
VARIANTS = {"crossfit_raw", "single_ibd", "crossfit_ibd"}


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


def bootstrap_single_split(
    edge_sets,
    *,
    train,
    evaluation,
    bandwidth: float,
    n_bootstrap: int,
    seed: int,
) -> SpeciesBootstrapResult:
    edge_map = {edges.species: edges for edges in edge_sets}
    prepared = prepare_transfer(
        [edge_map[name] for name in train],
        [edge_map[name] for name in evaluation],
        bandwidth=bandwidth,
    )
    observed = prepared.score(
        {name: edge_map[name].turnover for name in train},
        {name: edge_map[name].turnover for name in evaluation},
    )
    scores = np.asarray(
        [observed.species_scores[name] for name in evaluation if np.isfinite(observed.species_scores[name])],
        dtype=float,
    )
    bootstrap = centered_species_bootstrap_mean_test(scores, n_bootstrap=n_bootstrap, seed=seed)
    return SpeciesBootstrapResult(
        observed=TransferResult(
            statistic=float(scores.mean()),
            species_scores=observed.species_scores,
            n_eval_species=len(scores),
        ),
        species_scores=scores,
        null_statistics=bootstrap.null_means,
        null_studentized=bootstrap.null_studentized,
        p_value=bootstrap.p_value,
        null_mean=float(bootstrap.null_means.mean()),
        null_sd=float(bootstrap.null_means.std(ddof=1)),
        observed_studentized=bootstrap.observed_studentized,
    )


def main() -> int:
    parser = argparse.ArgumentParser(description="Compare small-panel TTF inference variants on frozen Iranian geometry.")
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
    parser.add_argument("--fold-seed", type=int, default=20260907)
    parser.add_argument("--ibd-degree", type=int, default=1)
    parser.add_argument("--seed", type=int, default=20260907)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    geometry = load_geometry(args.geometry)
    labels = tuple(sorted(geometry))
    train, evaluation = split_species(labels, eval_fraction=6.0 / len(labels), seed=args.split_seed)
    folds = balanced_species_folds(labels, n_folds=3, seed=args.fold_seed)

    p_values = []
    statistics = []
    null_means = []
    bandwidths = []
    for replicate in range(args.replicates):
        # World seed deliberately excludes the method variant: all variants see
        # the exact same synthetic worlds for paired method development.
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
        if args.variant in {"single_ibd", "crossfit_ibd"}:
            edge_sets = residualize_edge_sets(edge_sets, degree=args.ibd_degree)

        test_seed = seed_for(
            args.seed, "iran-ablation-bootstrap", args.variant,
            args.shared_fraction, args.amplitude, replicate,
        )
        if args.variant.startswith("crossfit"):
            result, used_folds = crossfit_species_bootstrap_test(
                edge_sets,
                n_folds=3,
                fold_seed=args.fold_seed,
                bandwidth=world.bandwidth,
                n_bootstrap=args.bootstrap,
                seed=test_seed,
            )
            if used_folds != folds:
                raise RuntimeError("cross-fit fold assignment drifted")
        else:
            result = bootstrap_single_split(
                edge_sets,
                train=train,
                evaluation=evaluation,
                bandwidth=world.bandwidth,
                n_bootstrap=args.bootstrap,
                seed=test_seed,
            )
        p_values.append(result.p_value)
        statistics.append(result.observed.statistic)
        null_means.append(result.null_mean)
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
        mean_statistic=float(np.mean(statistics)),
        mean_null_statistic=float(np.mean(null_means)),
        median_p_value=float(np.median(p)),
    )
    payload = {
        "schema": "ttf_iran_method_ablation_cell_v0.1",
        "variant": args.variant,
        "cell": asdict(cell),
        "design": {
            "geometry_schema": "ttf_iran_plateau_geometry_v0.1",
            "bandwidth": float(bandwidths[0]),
            "bandwidth_rule": "pooled median species-local kNN edge length after normalization",
            "k": int(args.k),
            "noise_sd": float(args.noise_sd),
            "transition_width": float(args.transition_width),
            "min_shared_crossing_species": int(args.min_crossing_species),
            "bootstrap_resamples": int(args.bootstrap),
            "single_train_species": list(train),
            "single_eval_species": list(evaluation),
            "crossfit_folds": [list(fold) for fold in folds],
            "ibd_residualization": "leave-one-edge-out rank-linear edge-length nuisance" if "ibd" in args.variant else None,
            "ibd_degree": int(args.ibd_degree) if "ibd" in args.variant else None,
            "master_seed": int(args.seed),
            "worlds_paired_across_variants": True,
            "empirical_genetic_outcomes_used": False,
            "known_boundary_labels_used": False,
        },
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"variant": args.variant, "cell": payload["cell"]}, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
