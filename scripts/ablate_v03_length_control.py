#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

from run_geometry_calibration import columns, load_geometry
from ttf.calibration import seed_for
from ttf.core import SpeciesSample, edge_turnover, spearman_rho, split_species
from ttf.geometry import simulate_fixed_geometry_boundary_world
from ttf.geometry_control import length_orthogonalized_turnover, partial_spearman_rho
from ttf.inference import centered_species_bootstrap_mean_test
from ttf.nulls import edges_on_fixed_graphs, fixed_graphs
from ttf.transfer import prepare_transfer

CANDIDATE_ORDER = (
    "eval_partial_only",
    "train_orthogonalized_only",
    "train_and_eval_conditioned",
)


def pack_train(prepared, values, *, residualize: bool, lengths):
    width = next(iter(values.values())).shape[1]
    n_edges = max(sl.stop for sl in prepared.train_slices.values())
    packed = np.empty((n_edges, width), dtype=float)
    for species in prepared.train_species:
        sl = prepared.train_slices[species]
        x = np.asarray(values[species], dtype=float)
        if residualize:
            for j in range(width):
                packed[sl, j] = length_orthogonalized_turnover(x[:, j], lengths[species])
        else:
            packed[sl, :] = x
    return packed


def score_variant(prepared, packed_train, eval_values, eval_lengths, *, partial: bool):
    width = packed_train.shape[1]
    scores: dict[str, np.ndarray] = {}
    stats = np.zeros(width, dtype=float)
    counts = np.zeros(width, dtype=int)
    for species in prepared.eval_species:
        pred = prepared.eval_projection[species] @ packed_train + prepared.eval_prior_offset[species][:, None]
        target = np.asarray(eval_values[species], dtype=float)
        length = np.asarray(eval_lengths[species], dtype=float)
        s = np.empty(width, dtype=float)
        for j in range(width):
            if partial:
                s[j] = partial_spearman_rho(pred[:, j], target[:, j], length)
            else:
                s[j] = spearman_rho(pred[:, j], target[:, j])
        scores[species] = s
        finite = np.isfinite(s)
        stats[finite] += s[finite]
        counts[finite] += 1
    if np.any(counts == 0):
        raise RuntimeError("one or more worlds have no finite held-out score")
    return stats / counts, scores


def p_for_column(scores, names, column, *, resamples, seed):
    x = np.asarray([scores[name][column] for name in names], dtype=float)
    x = x[np.isfinite(x)]
    return centered_species_bootstrap_mean_test(x, n_bootstrap=resamples, seed=seed).p_value


def main() -> int:
    p = argparse.ArgumentParser(description="Development-only v0.3 edge-length control ablation")
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
        raise RuntimeError("unexpected excluded species")
    labels = [g.species for g in geometries]
    split_seed = seed_for(args.seed, "fixed_geometry_split")
    train, evaluation = split_species(labels, eval_fraction=0.5, seed=split_seed)
    gm = {g.species: g for g in geometries}
    used = train + evaluation
    templates = [
        SpeciesSample(
            species=name,
            coordinates=gm[name].coordinates,
            trait=np.arange(len(gm[name].coordinates), dtype=float),
        )
        for name in used
    ]
    graphs = fixed_graphs(templates, k=args.k)
    edges = edges_on_fixed_graphs(templates, graphs)
    train_edges = [edges[name] for name in train]
    eval_edges = [edges[name] for name in evaluation]
    train_length = {name: edges[name].length for name in train}
    eval_length = {name: edges[name].length for name in evaluation}

    raw_prepared = prepare_transfer(
        train_edges, eval_edges,
        bandwidth=args.bandwidth,
        prior_strength=0.25,
        prior_mean=0.5,
        segment_points=5,
    )
    residual_prepared = prepare_transfer(
        train_edges, eval_edges,
        bandwidth=args.bandwidth,
        prior_strength=0.25,
        prior_mean=0.0,
        segment_points=5,
    )

    cells = [(0.0, 2.0), (0.0, 3.0), (1.0, 2.0)]
    result_by_variant = {
        name: {cell: {"p": [], "stat": []} for cell in cells}
        for name in ("original", *CANDIDATE_ORDER)
    }

    for shared, amplitude in cells:
        for batch_start in range(0, args.replicates, args.batch_size):
            stop = min(batch_start + args.batch_size, args.replicates)
            train_columns = {name: [] for name in train}
            eval_columns = {name: [] for name in evaluation}
            seeds = []
            for replicate in range(batch_start, stop):
                world = simulate_fixed_geometry_boundary_world(
                    geometries,
                    shared_fraction=shared,
                    amplitude=amplitude,
                    noise_sd=0.8,
                    transition_width=0.2,
                    seed=seed_for(args.seed, "fixed_geometry_world", shared, amplitude, replicate),
                )
                sm = {sample.species: sample for sample in world.samples}
                turnover = {name: edge_turnover(sm[name], graphs[name]) for name in used}
                for name in train:
                    train_columns[name].append(turnover[name])
                for name in evaluation:
                    eval_columns[name].append(turnover[name])
                seeds.append(seed_for(args.seed, "fixed_geometry_bootstrap", shared, amplitude, replicate))

            train_batch = {name: np.column_stack(train_columns[name]) for name in train}
            eval_batch = {name: np.column_stack(eval_columns[name]) for name in evaluation}
            raw_train = pack_train(raw_prepared, train_batch, residualize=False, lengths=train_length)
            residual_train = pack_train(residual_prepared, train_batch, residualize=True, lengths=train_length)

            variants = {
                "original": score_variant(raw_prepared, raw_train, eval_batch, eval_length, partial=False),
                "eval_partial_only": score_variant(raw_prepared, raw_train, eval_batch, eval_length, partial=True),
                "train_orthogonalized_only": score_variant(residual_prepared, residual_train, eval_batch, eval_length, partial=False),
                "train_and_eval_conditioned": score_variant(residual_prepared, residual_train, eval_batch, eval_length, partial=True),
            }
            for variant, (stats, score_map) in variants.items():
                store = result_by_variant[variant][(shared, amplitude)]
                for column, base_seed in enumerate(seeds):
                    pval = p_for_column(
                        score_map,
                        evaluation,
                        column,
                        resamples=args.resamples,
                        seed=seed_for(base_seed, "v03_ablation", variant),
                    )
                    store["p"].append(pval)
                    store["stat"].append(float(stats[column]))

    summaries = {}
    for variant, data in result_by_variant.items():
        cells_out = []
        for shared, amplitude in cells:
            pvals = np.asarray(data[(shared, amplitude)]["p"], dtype=float)
            stats = np.asarray(data[(shared, amplitude)]["stat"], dtype=float)
            cells_out.append({
                "shared_fraction": shared,
                "amplitude": amplitude,
                "rejection_rate": float(np.mean(pvals <= 0.05)),
                "mean_statistic": float(stats.mean()),
                "median_p": float(np.median(pvals)),
            })
        private = [c["rejection_rate"] for c in cells_out if c["shared_fraction"] == 0.0]
        power = [c["rejection_rate"] for c in cells_out if c["shared_fraction"] == 1.0][0]
        summaries[variant] = {
            "cells": cells_out,
            "development_max_private_rejection": float(max(private)),
            "development_power": float(power),
            "development_pass": bool(max(private) <= 0.10 and power >= 0.80),
        }

    selected = None
    for candidate in CANDIDATE_ORDER:
        if summaries[candidate]["development_pass"]:
            selected = candidate
            break

    payload = {
        "schema": "ttf_v03_length_control_ablation_v0.1",
        "status": "development_only_on_failed_gate_i_b_geometry",
        "claim_ready": False,
        "selection_rule_frozen_in_code_before_ablation_outcome": {
            "private_cells": [[0.0, 2.0], [0.0, 3.0]],
            "positive_cell": [1.0, 2.0],
            "max_private_rejection_ceiling": 0.10,
            "power_floor": 0.80,
            "candidate_preference_order": list(CANDIDATE_ORDER),
            "interpretation": "select the least invasive listed correction that satisfies both development point criteria",
        },
        "selected_candidate": selected,
        "fresh_external_qualification_required": True,
        "failed_panel_cannot_qualify_selected_candidate": True,
        "config": {
            "replicates": args.replicates,
            "resamples": args.resamples,
            "k": args.k,
            "bandwidth": args.bandwidth,
            "batch_size": args.batch_size,
            "seed": args.seed,
            "train_species": len(train),
            "eval_species": len(evaluation),
        },
        "variants": summaries,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"selected_candidate": selected, "variants": summaries}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
