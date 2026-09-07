#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

from run_geometry_calibration import columns, load_geometry
from ttf.batch import score_prepared_batch
from ttf.calibration import seed_for
from ttf.core import SpeciesSample, edge_turnover, spearman_rho, split_species
from ttf.geometry import simulate_fixed_geometry_boundary_world
from ttf.geometry_conditioned_batch import score_prepared_length_conditioned_batch
from ttf.inference import centered_species_bootstrap_mean_test
from ttf.nulls import edges_on_fixed_graphs, fixed_graphs
from ttf.transfer import prepare_transfer


def bootstrap_p(species_scores: np.ndarray, *, n_bootstrap: int, seed: int) -> float:
    finite = np.asarray(species_scores, dtype=float)
    finite = finite[np.isfinite(finite)]
    if len(finite) < 6:
        raise RuntimeError("fewer than six finite held-out species scores")
    return centered_species_bootstrap_mean_test(
        finite,
        n_bootstrap=n_bootstrap,
        seed=seed,
    ).p_value


def prediction_length_rho(
    prepared,
    train_turnover,
    eval_edge_length,
) -> dict[str, np.ndarray]:
    width = next(iter(train_turnover.values())).shape[1]
    n_train_edges = max(sl.stop for sl in prepared.train_slices.values())
    train_values = np.empty((n_train_edges, width), dtype=float)
    for species in prepared.train_species:
        sl = prepared.train_slices[species]
        train_values[sl, :] = np.asarray(train_turnover[species], dtype=float)

    out: dict[str, np.ndarray] = {}
    for species in prepared.eval_species:
        predicted = (
            prepared.eval_projection[species] @ train_values
            + prepared.eval_prior_offset[species][:, None]
        )
        length = np.asarray(eval_edge_length[species], dtype=float)
        out[species] = np.asarray(
            [spearman_rho(predicted[:, i], length) for i in range(width)],
            dtype=float,
        )
    return out


def mean_map_column(mapping: dict[str, np.ndarray], column: int) -> float:
    values = np.asarray(
        [mapping[name][column] for name in sorted(mapping)],
        dtype=float,
    )
    values = values[np.isfinite(values)]
    return float(values.mean()) if len(values) else float("nan")


def main() -> int:
    p = argparse.ArgumentParser(
        description=(
            "Development-only paired diagnostic of edge-length nuisance in the "
            "failed Gate-I-B geometry. This script cannot qualify v0.3."
        )
    )
    p.add_argument("--input", type=Path, required=True)
    p.add_argument("--coordinate-columns", type=columns, default=["x_km", "y_km", "z_km"])
    p.add_argument("--min-records", type=int, default=20)
    p.add_argument("--replicates", type=int, default=500)
    p.add_argument("--resamples", type=int, default=999)
    p.add_argument("--k", type=int, default=3)
    p.add_argument("--bandwidth", type=float, default=500.0)
    p.add_argument("--batch-size", type=int, default=25)
    p.add_argument("--seed", type=int, default=20260907)
    p.add_argument("--output", type=Path, required=True)
    args = p.parse_args()

    geometries, excluded = load_geometry(
        args.input,
        species_column="species",
        coordinate_columns=args.coordinate_columns,
        block_column=None,
        min_records=args.min_records,
    )
    if excluded:
        raise RuntimeError("diagnostic geometry unexpectedly excludes species")
    labels = [g.species for g in geometries]
    split_seed = seed_for(args.seed, "fixed_geometry_split")
    train, evaluation = split_species(labels, eval_fraction=0.5, seed=split_seed)

    geometry_map = {g.species: g for g in geometries}
    used = train + evaluation
    template_samples = [
        SpeciesSample(
            species=name,
            coordinates=geometry_map[name].coordinates,
            trait=np.arange(len(geometry_map[name].coordinates), dtype=float),
        )
        for name in used
    ]
    graphs = fixed_graphs(template_samples, k=args.k)
    template_edges = edges_on_fixed_graphs(template_samples, graphs)
    train_edges = [template_edges[name] for name in train]
    eval_edges = [template_edges[name] for name in evaluation]
    train_length = {name: template_edges[name].length for name in train}
    eval_length = {name: template_edges[name].length for name in evaluation}

    prepared_original = prepare_transfer(
        train_edges,
        eval_edges,
        bandwidth=args.bandwidth,
        prior_strength=0.25,
        prior_mean=0.5,
        segment_points=5,
    )
    prepared_conditioned = prepare_transfer(
        train_edges,
        eval_edges,
        bandwidth=args.bandwidth,
        prior_strength=0.25,
        prior_mean=0.0,
        segment_points=5,
    )

    cells = [(0.0, 2.0), (0.0, 3.0), (1.0, 2.0)]
    summaries: list[dict[str, object]] = []

    for shared, amplitude in cells:
        original_p: list[float] = []
        conditioned_p: list[float] = []
        original_stat: list[float] = []
        conditioned_stat: list[float] = []
        target_length_macro: list[float] = []
        original_prediction_length_macro: list[float] = []
        conditioned_prediction_length_macro: list[float] = []

        for batch_start in range(0, args.replicates, args.batch_size):
            batch_stop = min(batch_start + args.batch_size, args.replicates)
            train_columns = {name: [] for name in train}
            eval_columns = {name: [] for name in evaluation}
            null_seeds: list[int] = []

            for replicate in range(batch_start, batch_stop):
                world_seed = seed_for(
                    args.seed,
                    "fixed_geometry_world",
                    shared,
                    amplitude,
                    replicate,
                )
                null_seed = seed_for(
                    args.seed,
                    "fixed_geometry_bootstrap",
                    shared,
                    amplitude,
                    replicate,
                )
                world = simulate_fixed_geometry_boundary_world(
                    geometries,
                    shared_fraction=shared,
                    amplitude=amplitude,
                    noise_sd=0.8,
                    transition_width=0.2,
                    seed=world_seed,
                )
                sample_map = {sample.species: sample for sample in world.samples}
                turnover = {
                    name: edge_turnover(sample_map[name], graphs[name]) for name in used
                }
                for name in train:
                    train_columns[name].append(turnover[name])
                for name in evaluation:
                    eval_columns[name].append(turnover[name])
                null_seeds.append(null_seed)

            train_batch = {
                name: np.column_stack(train_columns[name]) for name in train
            }
            eval_batch = {
                name: np.column_stack(eval_columns[name]) for name in evaluation
            }
            original = score_prepared_batch(
                prepared_original,
                train_batch,
                eval_batch,
            )
            conditioned = score_prepared_length_conditioned_batch(
                prepared_conditioned,
                train_batch,
                eval_batch,
                train_edge_length=train_length,
                eval_edge_length=eval_length,
            )
            original_prediction_length = prediction_length_rho(
                prepared_original,
                train_batch,
                eval_length,
            )

            for column, null_seed in enumerate(null_seeds):
                raw_species = np.asarray(
                    [original.species_scores[name][column] for name in evaluation],
                    dtype=float,
                )
                conditioned_species = np.asarray(
                    [conditioned.species_scores[name][column] for name in evaluation],
                    dtype=float,
                )
                original_p.append(
                    bootstrap_p(raw_species, n_bootstrap=args.resamples, seed=null_seed)
                )
                conditioned_p.append(
                    bootstrap_p(
                        conditioned_species,
                        n_bootstrap=args.resamples,
                        seed=seed_for(null_seed, "length_conditioned"),
                    )
                )
                original_stat.append(float(original.statistics[column]))
                conditioned_stat.append(float(conditioned.statistics[column]))
                target_length_macro.append(
                    mean_map_column(conditioned.eval_turnover_length_rho, column)
                )
                original_prediction_length_macro.append(
                    mean_map_column(original_prediction_length, column)
                )
                conditioned_prediction_length_macro.append(
                    mean_map_column(conditioned.eval_prediction_length_rho, column)
                )

        raw_p = np.asarray(original_p, dtype=float)
        adj_p = np.asarray(conditioned_p, dtype=float)
        summaries.append(
            {
                "shared_fraction": shared,
                "amplitude": amplitude,
                "worlds": args.replicates,
                "original": {
                    "rejection_rate": float(np.mean(raw_p <= 0.05)),
                    "mean_statistic": float(np.mean(original_stat)),
                    "median_p": float(np.median(raw_p)),
                },
                "length_conditioned_v03_candidate": {
                    "rejection_rate": float(np.mean(adj_p <= 0.05)),
                    "mean_statistic": float(np.mean(conditioned_stat)),
                    "median_p": float(np.median(adj_p)),
                },
                "geometry_diagnostic": {
                    "mean_eval_turnover_vs_edge_length_spearman": float(np.mean(target_length_macro)),
                    "mean_original_prediction_vs_edge_length_spearman": float(np.mean(original_prediction_length_macro)),
                    "mean_conditioned_prediction_vs_edge_length_spearman": float(np.mean(conditioned_prediction_length_macro)),
                },
            }
        )

    payload = {
        "schema": "ttf_v03_gate_i_b_length_nuisance_diagnostic_v0.1",
        "status": "development_only_on_failed_gate_i_b_geometry",
        "claim_ready": False,
        "may_select_successor_candidate": True,
        "may_qualify_successor": False,
        "failed_geometry_reuse_rule": (
            "The failed Gate-I-B panel may diagnose and develop a successor but may "
            "not provide confirmatory qualification for that successor. Any selected "
            "v0.3 statistic requires prospectively frozen validation on fresh geometry."
        ),
        "geometry_input": str(args.input),
        "split_seed": split_seed,
        "train_species": len(train),
        "eval_species": len(evaluation),
        "config": {
            "k": args.k,
            "bandwidth": args.bandwidth,
            "worlds_per_cell": args.replicates,
            "bootstrap_resamples": args.resamples,
            "batch_size": args.batch_size,
            "seed": args.seed,
            "candidate": "train_length_orthogonalized_turnover_plus_eval_partial_spearman_given_edge_length_rank",
            "candidate_prior_mean": 0.0,
            "candidate_prior_strength": 0.25,
        },
        "cells": summaries,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    print(json.dumps(summaries, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
