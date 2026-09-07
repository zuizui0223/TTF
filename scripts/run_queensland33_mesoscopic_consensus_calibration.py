#!/usr/bin/env python3
from __future__ import annotations

import argparse
from dataclasses import asdict
import json
from pathlib import Path

import numpy as np

from run_queensland33_mesoscopic_calibration import (
    candidate_cuts,
    geometry_graphs,
    load_mesoscopic_layout,
)
from run_queensland33_topology_calibration import simulate_ordered_world
from ttf.calibration import CalibrationCell, seed_for
from ttf.core import split_species
from ttf.mesoscopic_consensus import consensus_bootstrap_test
from ttf.nulls import edges_on_fixed_graphs


def main() -> int:
    parser = argparse.ArgumentParser(
        description=(
            "Calibrate the prospective species-normalized mesoscopic consensus "
            "TTF estimator on the frozen Queensland topology envelope."
        )
    )
    parser.add_argument("--layouts", type=Path, required=True)
    parser.add_argument("--layout-index", type=int, required=True)
    parser.add_argument("--shared-fraction", type=float, required=True)
    parser.add_argument("--amplitude", type=float, required=True)
    parser.add_argument("--replicates", type=int, default=100)
    parser.add_argument("--bootstrap", type=int, default=999)
    parser.add_argument("--alpha", type=float, default=0.05)
    parser.add_argument("--k", type=int, default=2)
    parser.add_argument("--noise-sd", type=float, default=0.8)
    parser.add_argument("--transition-width", type=float, default=0.20)
    parser.add_argument("--min-shared-crossing-species", type=int, default=14)
    parser.add_argument("--prior-strength", type=float, default=1.0)
    parser.add_argument("--evidence-temperature", type=float, default=1.0)
    parser.add_argument("--split-seed", type=int, default=20260907)
    parser.add_argument("--seed", type=int, default=20260911)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    if not np.isclose(float(args.evidence_temperature), 1.0, atol=0.0, rtol=0.0):
        raise ValueError("Queensland v0.3b qualification prospectively fixes evidence_temperature=1.0")

    geometry, layout_meta = load_mesoscopic_layout(args.layouts, int(args.layout_index))
    labels = tuple(sorted(geometry))
    train, evaluation = split_species(
        labels,
        eval_fraction=11.0 / len(labels),
        seed=int(args.split_seed),
    )
    if len(train) != 10 or len(evaluation) != 11:
        raise RuntimeError("Queensland consensus benchmark requires frozen 10/11 split")

    cuts = candidate_cuts(geometry)
    if len(cuts) != 7:
        raise RuntimeError("Queensland consensus benchmark requires seven inter-basin cuts")
    graphs = geometry_graphs(geometry, k=int(args.k))

    p_values: list[float] = []
    statistics: list[float] = []
    max_field_probability: list[float] = []
    field_entropy: list[float] = []
    finite_eval: list[int] = []
    informative_train: list[int] = []
    optimizer_iterations: list[int] = []

    for replicate in range(int(args.replicates)):
        samples, _bandwidth, _shared_cut, _crossing = simulate_ordered_world(
            geometry,
            shared_fraction=float(args.shared_fraction),
            amplitude=float(args.amplitude),
            noise_sd=float(args.noise_sd),
            transition_width=float(args.transition_width),
            min_shared_crossing_species=int(args.min_shared_crossing_species),
            k=int(args.k),
            seed=seed_for(
                int(args.seed),
                "qld33-consensus-world",
                int(args.layout_index),
                float(args.shared_fraction),
                float(args.amplitude),
                replicate,
            ),
        )
        edge_map = edges_on_fixed_graphs(samples, graphs)
        result = consensus_bootstrap_test(
            [edge_map[s] for s in train],
            [edge_map[s] for s in evaluation],
            cuts,
            prior_strength=float(args.prior_strength),
            evidence_temperature=float(args.evidence_temperature),
            n_bootstrap=int(args.bootstrap),
            seed=seed_for(
                int(args.seed),
                "qld33-consensus-bootstrap",
                int(args.layout_index),
                float(args.shared_fraction),
                float(args.amplitude),
                replicate,
            ),
        )
        p_values.append(float(result.p_value))
        statistics.append(float(result.observed.statistic))
        probability = np.asarray(result.observed.field.probability, dtype=float)
        max_field_probability.append(float(np.max(probability)))
        field_entropy.append(float(-np.sum(probability * np.log(probability))))
        finite_eval.append(int(result.observed.n_eval_species))
        informative_train.append(int(result.observed.field.n_informative_species))
        optimizer_iterations.append(int(result.observed.field.iterations))

    p = np.asarray(p_values, dtype=float)
    cell = CalibrationCell(
        shared_fraction=float(args.shared_fraction),
        amplitude=float(args.amplitude),
        n_replicates=int(args.replicates),
        alpha=float(args.alpha),
        rejection_rate=float(np.mean(p <= float(args.alpha))),
        mean_statistic=float(np.mean(statistics)),
        mean_null_statistic=0.0,
        median_p_value=float(np.median(p)),
    )
    payload = {
        "schema": "ttf_queensland33_mesoscopic_consensus_cell_v0.3b",
        "layout_index": int(args.layout_index),
        "cell": asdict(cell),
        "design": {
            "estimator": "species-normalized cut posterior plus observability-conditioned consensus held-out log-score gain",
            "n_candidate_cuts": int(len(cuts)),
            "n_eligible_species": int(len(labels)),
            "train_species": list(train),
            "eval_species": list(evaluation),
            "split": "10 train / 11 evaluation species",
            "split_seed": int(args.split_seed),
            "layout_schema": layout_meta["layout_schema"],
            "layout_mode": layout_meta["layout_mode"],
            "exact_frozen_layout_envelope": True,
            "layout_selected_by_outcome": False,
            "named_Mary_Brisbane_boundary_used": False,
            "genetic_outcomes_used": False,
            "k": int(args.k),
            "noise_sd": float(args.noise_sd),
            "transition_width": float(args.transition_width),
            "min_shared_crossing_species": int(args.min_shared_crossing_species),
            "prior_strength": float(args.prior_strength),
            "evidence_temperature": float(args.evidence_temperature),
            "bootstrap_resamples": int(args.bootstrap),
            "master_seed": int(args.seed),
            "candidate_registered_before_v0.3a_result": True,
        },
        "diagnostics": {
            "mean_max_field_probability": float(np.mean(max_field_probability)),
            "mean_field_entropy": float(np.mean(field_entropy)),
            "mean_finite_eval_species": float(np.mean(finite_eval)),
            "min_finite_eval_species": int(min(finite_eval)),
            "mean_informative_train_species": float(np.mean(informative_train)),
            "min_informative_train_species": int(min(informative_train)),
            "mean_optimizer_iterations": float(np.mean(optimizer_iterations)),
            "max_optimizer_iterations": int(max(optimizer_iterations)),
        },
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    print(
        json.dumps(
            {
                "layout": int(args.layout_index),
                "shared_fraction": float(args.shared_fraction),
                "amplitude": float(args.amplitude),
                "rejection_rate": cell.rejection_rate,
                "mean_statistic": cell.mean_statistic,
                "mean_max_field_probability": payload["diagnostics"]["mean_max_field_probability"],
            },
            sort_keys=True,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
