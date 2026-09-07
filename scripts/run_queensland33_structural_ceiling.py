#!/usr/bin/env python3
from __future__ import annotations

import argparse
from dataclasses import asdict
import json
from pathlib import Path

import numpy as np

from ttf.calibration import CalibrationCell, seed_for
from ttf.core import split_species
from ttf.inference import heldout_species_bootstrap_test
from run_queensland33_detection_floor import load_frozen_layout
from run_queensland33_topology_calibration import simulate_ordered_world


def main() -> int:
    parser = argparse.ArgumentParser(description="Estimate the noise-free structural identifiability ceiling on frozen Queensland layouts.")
    parser.add_argument("--layouts", type=Path, required=True)
    parser.add_argument("--layout-index", type=int, required=True)
    parser.add_argument("--shared-fraction", type=float, choices=(0.0, 1.0), required=True)
    parser.add_argument("--replicates", type=int, default=100)
    parser.add_argument("--bootstrap", type=int, default=1999)
    parser.add_argument("--alpha", type=float, default=0.05)
    parser.add_argument("--k", type=int, default=2)
    parser.add_argument("--transition-width", type=float, default=0.20)
    parser.add_argument("--min-shared-crossing-species", type=int, default=14)
    parser.add_argument("--split-seed", type=int, default=20260907)
    parser.add_argument("--seed", type=int, default=20260909)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    geometry = load_frozen_layout(args.layouts, int(args.layout_index))
    labels = tuple(sorted(geometry))
    train, evaluation = split_species(labels, eval_fraction=11.0 / len(labels), seed=int(args.split_seed))
    if len(train) != 10 or len(evaluation) != 11:
        raise RuntimeError("Queensland ceiling requires frozen 10/11 split")

    p_values: list[float] = []
    statistics: list[float] = []
    null_means: list[float] = []
    train_crossing: list[int] = []
    eval_crossing: list[int] = []
    total_crossing: list[int] = []

    for replicate in range(int(args.replicates)):
        samples, bandwidth, _shared_cut, crossing = simulate_ordered_world(
            geometry,
            shared_fraction=float(args.shared_fraction),
            amplitude=1.0,
            noise_sd=0.0,
            transition_width=float(args.transition_width),
            min_shared_crossing_species=int(args.min_shared_crossing_species),
            k=int(args.k),
            seed=seed_for(args.seed, "qld33-structural-ceiling-world", args.layout_index, args.shared_fraction, replicate),
        )
        result = heldout_species_bootstrap_test(
            samples,
            train_species=train,
            eval_species=evaluation,
            k=int(args.k),
            bandwidth=float(bandwidth),
            n_bootstrap=int(args.bootstrap),
            seed=seed_for(args.seed, "qld33-structural-ceiling-bootstrap", args.layout_index, args.shared_fraction, replicate),
        )
        p_values.append(float(result.p_value))
        statistics.append(float(result.observed.statistic))
        null_means.append(float(result.null_mean))
        if float(args.shared_fraction) == 1.0:
            crossing_set = set(crossing)
            train_crossing.append(sum(name in crossing_set for name in train))
            eval_crossing.append(sum(name in crossing_set for name in evaluation))
            total_crossing.append(len(crossing_set))

    p = np.asarray(p_values, dtype=float)
    cell = CalibrationCell(
        shared_fraction=float(args.shared_fraction),
        amplitude=1.0,
        n_replicates=int(args.replicates),
        alpha=float(args.alpha),
        rejection_rate=float(np.mean(p <= float(args.alpha))),
        mean_statistic=float(np.mean(statistics)),
        mean_null_statistic=float(np.mean(null_means)),
        median_p_value=float(np.median(p)),
    )

    coverage = None
    if eval_crossing:
        coverage = {
            "mean_total_crossing_species": float(np.mean(total_crossing)),
            "min_total_crossing_species": int(min(total_crossing)),
            "mean_train_crossing_species": float(np.mean(train_crossing)),
            "min_train_crossing_species": int(min(train_crossing)),
            "mean_eval_crossing_species": float(np.mean(eval_crossing)),
            "min_eval_crossing_species": int(min(eval_crossing)),
            "mean_eval_crossing_fraction": float(np.mean(np.asarray(eval_crossing) / len(evaluation))),
        }

    payload = {
        "schema": "ttf_queensland33_structural_ceiling_cell_v0.1",
        "layout_index": int(args.layout_index),
        "cell": asdict(cell),
        "coverage": coverage,
        "design": {
            "frozen_layout_schema": "ttf_queensland33_topology_layouts_frozen_v0.2",
            "n_eligible_species": len(labels),
            "split": "10 train / 11 evaluation species",
            "train_species": list(train),
            "eval_species": list(evaluation),
            "split_seed": int(args.split_seed),
            "k": int(args.k),
            "amplitude": 1.0,
            "noise_sd": 0.0,
            "interpretation_of_amplitude": "irrelevant scale under within-species ranks when noise=0; this is the practical infinite-SNR ceiling",
            "transition_width": float(args.transition_width),
            "min_shared_crossing_species": int(args.min_shared_crossing_species),
            "bootstrap_resamples": int(args.bootstrap),
            "master_seed": int(args.seed),
            "estimator_changed": False,
            "genetic_outcomes_used": False,
            "named_Mary_Brisbane_boundary_used": False,
        },
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    print(json.dumps({
        "layout": int(args.layout_index),
        "shared_fraction": float(args.shared_fraction),
        "rejection_rate": cell.rejection_rate,
        "mean_statistic": cell.mean_statistic,
        "coverage": coverage,
    }, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
