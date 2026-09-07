#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

from ttf.calibration import seed_for
from ttf.core import build_species_edges, spearman_rho, split_species
from ttf.inference import centered_species_bootstrap_mean_test
from ttf.semisynthetic import normalize_geometry
from run_queensland33_detection_floor import load_frozen_layout
from run_queensland33_topology_calibration import BASINS, crossing_species, simulate_ordered_world


def eligible_normalized_cuts(geometry: dict[str, np.ndarray], *, min_total_crossing: int = 14):
    normalized = normalize_geometry(geometry)
    raw_cuts = np.arange(0.5, len(BASINS) - 0.5, 1.0)
    cuts = (raw_cuts - float(normalized.center[0])) / float(normalized.scale)
    norm_geometry = {
        s: np.column_stack([
            (geometry[s][:, 0] - normalized.center[0]) / normalized.scale,
            np.zeros(len(geometry[s])),
        ])
        for s in geometry
    }
    out = []
    for cut in cuts:
        crossing = crossing_species(norm_geometry, float(cut))
        if len(crossing) >= int(min_total_crossing):
            out.append((float(cut), crossing))
    if not out:
        raise RuntimeError("no geometry-only oracle cut satisfies crossing requirement")
    return norm_geometry, out


def oracle_segment_exposure(start: np.ndarray, end: np.ndarray, *, cut: float, bandwidth: float, n_points: int = 9) -> np.ndarray:
    t = np.linspace(0.0, 1.0, int(n_points))
    points = start[:, None, :] + t[None, :, None] * (end - start)[:, None, :]
    dx = (points[:, :, 0] - float(cut)) / float(bandwidth)
    return np.exp(-0.5 * dx * dx).mean(axis=1)


def score_oracle(samples, *, eval_species, cut: float, bandwidth: float, k: int, crossing_only: bool, crossing_set: set[str]):
    by_species = {sample.species: sample for sample in samples}
    scores = []
    used = []
    for species in eval_species:
        if crossing_only and species not in crossing_set:
            continue
        edges = build_species_edges(by_species[species], k=int(k))
        exposure = oracle_segment_exposure(edges.start, edges.end, cut=float(cut), bandwidth=float(bandwidth))
        score = spearman_rho(exposure, edges.turnover)
        if np.isfinite(score):
            scores.append(float(score))
            used.append(species)
    return np.asarray(scores, dtype=float), tuple(used)


def main() -> int:
    parser = argparse.ArgumentParser(description="Separate field-learning from held-out scoring using a true-boundary oracle on Queensland layouts.")
    parser.add_argument("--layouts", type=Path, required=True)
    parser.add_argument("--layout-index", type=int, required=True)
    parser.add_argument("--shared-fraction", type=float, choices=(0.0, 1.0), required=True)
    parser.add_argument("--replicates", type=int, default=100)
    parser.add_argument("--bootstrap", type=int, default=1999)
    parser.add_argument("--alpha", type=float, default=0.05)
    parser.add_argument("--k", type=int, default=2)
    parser.add_argument("--transition-width", type=float, default=0.20)
    parser.add_argument("--split-seed", type=int, default=20260907)
    parser.add_argument("--seed", type=int, default=20260910)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    geometry = load_frozen_layout(args.layouts, int(args.layout_index))
    labels = tuple(sorted(geometry))
    train, evaluation = split_species(labels, eval_fraction=11.0 / len(labels), seed=int(args.split_seed))
    normalized_geometry, possible_cuts = eligible_normalized_cuts(geometry, min_total_crossing=14)

    all_p = []
    crossing_p = []
    all_mean = []
    crossing_mean = []
    crossing_n = []

    for replicate in range(int(args.replicates)):
        samples, bandwidth, shared_cut, shared_crossing = simulate_ordered_world(
            geometry,
            shared_fraction=float(args.shared_fraction),
            amplitude=1.0,
            noise_sd=0.0,
            transition_width=float(args.transition_width),
            min_shared_crossing_species=14,
            k=int(args.k),
            seed=seed_for(args.seed, "qld33-oracle-world", args.layout_index, args.shared_fraction, replicate),
        )
        if float(args.shared_fraction) == 1.0:
            if shared_cut is None:
                raise RuntimeError("shared oracle world is missing its true cut")
            oracle_cut = float(shared_cut)
            oracle_crossing = set(shared_crossing)
        else:
            # Geometry-only arbitrary cut, drawn independently of every private trait transition.
            rng = np.random.default_rng(seed_for(args.seed, "qld33-oracle-null-cut", args.layout_index, replicate))
            candidates = []
            for cut, crossing in possible_cuts:
                n_eval = sum(species in set(crossing) for species in evaluation)
                if n_eval >= 6:
                    candidates.append((cut, crossing))
            if not candidates:
                raise RuntimeError("no geometry-only null cut has >=6 evaluable held-out species")
            idx = int(rng.integers(0, len(candidates)))
            oracle_cut, crossing = candidates[idx]
            oracle_crossing = set(crossing)

        scores_all, _ = score_oracle(
            samples,
            eval_species=evaluation,
            cut=oracle_cut,
            bandwidth=float(bandwidth),
            k=int(args.k),
            crossing_only=False,
            crossing_set=oracle_crossing,
        )
        scores_cross, used_cross = score_oracle(
            samples,
            eval_species=evaluation,
            cut=oracle_cut,
            bandwidth=float(bandwidth),
            k=int(args.k),
            crossing_only=True,
            crossing_set=oracle_crossing,
        )
        if len(scores_all) < 6 or len(scores_cross) < 6:
            raise RuntimeError("oracle diagnostic has fewer than six finite species scores")
        test_all = centered_species_bootstrap_mean_test(
            scores_all,
            n_bootstrap=int(args.bootstrap),
            seed=seed_for(args.seed, "qld33-oracle-all-bootstrap", args.layout_index, args.shared_fraction, replicate),
        )
        test_cross = centered_species_bootstrap_mean_test(
            scores_cross,
            n_bootstrap=int(args.bootstrap),
            seed=seed_for(args.seed, "qld33-oracle-cross-bootstrap", args.layout_index, args.shared_fraction, replicate),
        )
        all_p.append(test_all.p_value)
        crossing_p.append(test_cross.p_value)
        all_mean.append(test_all.observed_mean)
        crossing_mean.append(test_cross.observed_mean)
        crossing_n.append(len(used_cross))

    all_p = np.asarray(all_p, dtype=float)
    crossing_p = np.asarray(crossing_p, dtype=float)
    payload = {
        "schema": "ttf_queensland33_oracle_diagnostic_cell_v0.1",
        "layout_index": int(args.layout_index),
        "shared_fraction": float(args.shared_fraction),
        "oracle_all": {
            "rejection_rate": float(np.mean(all_p <= float(args.alpha))),
            "mean_species_score": float(np.mean(all_mean)),
            "median_p_value": float(np.median(all_p)),
        },
        "oracle_crossing": {
            "rejection_rate": float(np.mean(crossing_p <= float(args.alpha))),
            "mean_species_score": float(np.mean(crossing_mean)),
            "median_p_value": float(np.median(crossing_p)),
            "mean_n_species": float(np.mean(crossing_n)),
            "min_n_species": int(min(crossing_n)),
        },
        "design": {
            "diagnostic_only": True,
            "true_cut_used_only_when_shared_fraction_is_1": True,
            "private_null_cut_independent_of_trait": True,
            "private_null_cut_geometry_only": True,
            "noise_sd": 0.0,
            "amplitude": 1.0,
            "transition_width": float(args.transition_width),
            "same_frozen_10_11_split": True,
            "k": int(args.k),
            "bootstrap_resamples": int(args.bootstrap),
            "genetic_outcomes_used": False,
            "named_Mary_Brisbane_boundary_used": False,
        },
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    print(json.dumps({
        "layout": int(args.layout_index),
        "shared_fraction": float(args.shared_fraction),
        "oracle_all": payload["oracle_all"],
        "oracle_crossing": payload["oracle_crossing"],
    }, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
