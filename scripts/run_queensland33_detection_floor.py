#!/usr/bin/env python3
from __future__ import annotations

import argparse
from dataclasses import asdict
import json
from pathlib import Path

import numpy as np

from ttf.calibration import CalibrationCell, seed_for
from ttf.core import SpeciesSample, split_species
from ttf.inference import heldout_species_bootstrap_test
from run_queensland33_topology_calibration import simulate_ordered_world

BASINS = ("GOL", "LOG", "BRI", "PIN", "MCY", "NOO", "TCB", "MRY")


def load_frozen_layout(path: Path, layout_index: int) -> dict[str, np.ndarray]:
    payload = json.loads(path.read_text())
    if payload.get("schema") != "ttf_queensland33_topology_layouts_frozen_v0.2":
        raise RuntimeError("unexpected frozen Queensland layout schema")
    if tuple(payload["basin_column_order"]) != BASINS:
        raise RuntimeError("basin order drift")
    layouts = {int(item["i"]): item for item in payload["layouts"]}
    if int(layout_index) not in layouts:
        raise RuntimeError(f"unknown layout {layout_index}")
    row = layouts[int(layout_index)]
    geometry: dict[str, np.ndarray] = {}
    for species, counts in row["c"].items():
        if len(counts) != len(BASINS):
            raise RuntimeError(f"basin vector drift for {species}")
        supported = [i for i, n in enumerate(counts) if int(n) >= 3]
        if len(supported) < 4:
            raise RuntimeError(f"eligibility drift for {species}")
        x = np.asarray(supported, dtype=float)
        geometry[str(species)] = np.column_stack([x, np.zeros_like(x)])
    if len(geometry) != 21:
        raise RuntimeError(f"expected 21 eligible species, got {len(geometry)}")
    return geometry


def main() -> int:
    parser = argparse.ArgumentParser(description="Estimate the Queensland topology detection floor on the exact frozen 12-layout envelope.")
    parser.add_argument("--layouts", type=Path, required=True)
    parser.add_argument("--layout-index", type=int, required=True)
    parser.add_argument("--shared-fraction", type=float, required=True)
    parser.add_argument("--amplitude", type=float, required=True)
    parser.add_argument("--replicates", type=int, default=50)
    parser.add_argument("--bootstrap", type=int, default=999)
    parser.add_argument("--alpha", type=float, default=0.05)
    parser.add_argument("--k", type=int, default=2)
    parser.add_argument("--noise-sd", type=float, default=0.8)
    parser.add_argument("--transition-width", type=float, default=0.20)
    parser.add_argument("--min-shared-crossing-species", type=int, default=14)
    parser.add_argument("--split-seed", type=int, default=20260907)
    parser.add_argument("--seed", type=int, default=20260907)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    if (float(args.shared_fraction), float(args.amplitude)) not in {
        (0.0, 4.0), (1.0, 2.5), (1.0, 3.0), (1.0, 4.0)
    }:
        raise ValueError("detection-floor runner accepts only the frozen scan cells")

    geometry = load_frozen_layout(args.layouts, args.layout_index)
    labels = tuple(sorted(geometry))
    train, evaluation = split_species(labels, eval_fraction=11.0 / len(labels), seed=args.split_seed)
    if len(train) != 10 or len(evaluation) != 11:
        raise RuntimeError("Queensland benchmark requires frozen 10/11 split")

    p_values = []
    statistics = []
    null_means = []
    bandwidths = []
    crossing_counts = []
    for replicate in range(int(args.replicates)):
        samples, bandwidth, _cut, crossing = simulate_ordered_world(
            geometry,
            shared_fraction=float(args.shared_fraction),
            amplitude=float(args.amplitude),
            noise_sd=float(args.noise_sd),
            transition_width=float(args.transition_width),
            min_shared_crossing_species=int(args.min_shared_crossing_species),
            k=int(args.k),
            seed=seed_for(args.seed, "qld33-detection-floor-world", args.layout_index, args.shared_fraction, args.amplitude, replicate),
        )
        result = heldout_species_bootstrap_test(
            samples,
            train_species=train,
            eval_species=evaluation,
            k=int(args.k),
            bandwidth=float(bandwidth),
            n_bootstrap=int(args.bootstrap),
            seed=seed_for(args.seed, "qld33-detection-floor-bootstrap", args.layout_index, args.shared_fraction, args.amplitude, replicate),
        )
        p_values.append(result.p_value)
        statistics.append(result.observed.statistic)
        null_means.append(result.null_mean)
        bandwidths.append(float(bandwidth))
        crossing_counts.append(len(crossing))

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
        "schema": "ttf_queensland33_detection_floor_cell_v0.1",
        "layout_index": int(args.layout_index),
        "cell": asdict(cell),
        "design": {
            "frozen_layout_schema": "ttf_queensland33_topology_layouts_frozen_v0.2",
            "layout_selected_by_outcome": False,
            "genetic_outcomes_used": False,
            "named_Mary_Brisbane_boundary_used": False,
            "n_eligible_species": len(labels),
            "train_species": list(train),
            "eval_species": list(evaluation),
            "split": "10 train / 11 evaluation species",
            "split_seed": int(args.split_seed),
            "k": int(args.k),
            "noise_sd": float(args.noise_sd),
            "transition_width": float(args.transition_width),
            "min_shared_crossing_species": int(args.min_shared_crossing_species),
            "bootstrap_resamples": int(args.bootstrap),
            "master_seed": int(args.seed),
            "mean_bandwidth": float(np.mean(bandwidths)),
            "world_seed_family": "qld33-detection-floor-world",
        },
        "diagnostics": {
            "mean_crossing_species": float(np.mean(crossing_counts)) if crossing_counts else 0.0,
            "min_crossing_species": int(min(crossing_counts)) if crossing_counts else 0,
            "max_crossing_species": int(max(crossing_counts)) if crossing_counts else 0,
        },
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"layout": args.layout_index, "cell": payload["cell"], "diagnostics": payload["diagnostics"]}, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
