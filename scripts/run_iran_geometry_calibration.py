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
from ttf.semisynthetic import simulate_on_geometry

EARTH_RADIUS_KM = 6371.0088


def load_geometry(path: Path) -> dict[str, np.ndarray]:
    payload = json.loads(path.read_text())
    if payload.get("schema") != "ttf_iran_plateau_geometry_v0.1":
        raise RuntimeError("unexpected Iranian geometry schema")
    geometry: dict[str, np.ndarray] = {}
    all_latlon: list[tuple[float, float]] = []
    raw: dict[str, np.ndarray] = {}
    for item in payload["species"]:
        species = str(item["species"])
        latlon = np.asarray(
            [[row["latitude"], row["longitude"]] for row in item["populations"]],
            dtype=float,
        )
        raw[species] = latlon
        all_latlon.extend(map(tuple, latlon.tolist()))

    pooled = np.asarray(all_latlon, dtype=float)
    lat0 = float(np.deg2rad(pooled[:, 0].mean()))
    lon0 = float(np.deg2rad(pooled[:, 1].mean()))
    lat_center = float(np.deg2rad(pooled[:, 0].mean()))
    for species, latlon in raw.items():
        lat = np.deg2rad(latlon[:, 0])
        lon = np.deg2rad(latlon[:, 1])
        x = EARTH_RADIUS_KM * np.cos(lat0) * (lon - lon0)
        y = EARTH_RADIUS_KM * (lat - lat_center)
        geometry[species] = np.column_stack([x, y])
    return geometry


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Run TTF sharedness calibration on frozen Iranian Plateau sampling geometry."
    )
    parser.add_argument("--geometry", type=Path, required=True)
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
    parser.add_argument("--seed", type=int, default=20260907)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    if args.replicates < 1 or args.bootstrap < 99:
        raise ValueError("replicates must be positive and bootstrap >= 99")
    geometry = load_geometry(args.geometry)
    labels = tuple(sorted(geometry))
    # Nine species cannot support the v0.2 50/50 split because inference requires
    # at least six held-out species.  Freeze a 3-train / 6-evaluation split using
    # only taxon names and a prospective seed.
    train, evaluation = split_species(
        labels,
        eval_fraction=6.0 / len(labels),
        seed=args.split_seed,
    )
    if len(train) != 3 or len(evaluation) != 6:
        raise RuntimeError("Iranian benchmark did not resolve to a 3/6 species split")

    p_values: list[float] = []
    statistics: list[float] = []
    null_means: list[float] = []
    bandwidths: list[float] = []
    crossing_counts: list[int] = []
    for replicate in range(args.replicates):
        world = simulate_on_geometry(
            geometry,
            shared_fraction=args.shared_fraction,
            amplitude=args.amplitude,
            noise_sd=args.noise_sd,
            transition_width=args.transition_width,
            min_shared_crossing_species=args.min_crossing_species,
            k=args.k,
            seed=seed_for(args.seed, "iran-geometry-world", args.shared_fraction, args.amplitude, replicate),
        )
        result = heldout_species_bootstrap_test(
            world.samples,
            train_species=train,
            eval_species=evaluation,
            k=args.k,
            bandwidth=world.bandwidth,
            n_bootstrap=args.bootstrap,
            seed=seed_for(args.seed, "iran-geometry-bootstrap", args.shared_fraction, args.amplitude, replicate),
        )
        p_values.append(result.p_value)
        statistics.append(result.observed.statistic)
        null_means.append(result.null_mean)
        bandwidths.append(world.bandwidth)
        crossing_counts.append(len(world.crossing_species))

    if not np.allclose(bandwidths, bandwidths[0], atol=1e-12, rtol=0.0):
        raise RuntimeError("outcome-free geometry bandwidth drifted across worlds")
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
        "schema": "ttf_iran_geometry_calibration_cell_v0.1",
        "cell": asdict(cell),
        "design": {
            "inference": "heldout_species_bootstrap",
            "geometry_schema": "ttf_iran_plateau_geometry_v0.1",
            "coordinate_projection": "local equirectangular km; pooled RMS normalization inside TTF",
            "bandwidth_rule": "pooled median species-local kNN edge length after normalization",
            "bandwidth": float(bandwidths[0]),
            "k": int(args.k),
            "noise_sd": float(args.noise_sd),
            "transition_width": float(args.transition_width),
            "min_shared_crossing_species": int(args.min_crossing_species),
            "split_seed": int(args.split_seed),
            "master_seed": int(args.seed),
            "train_species": list(train),
            "eval_species": list(evaluation),
            "bootstrap_resamples": int(args.bootstrap),
            "outcome_values_used": False,
            "known_boundary_labels_used": False,
        },
        "crossing_species_count": {
            "min": int(min(crossing_counts)),
            "median": float(np.median(crossing_counts)),
            "max": int(max(crossing_counts)),
        },
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"cell": payload["cell"], "bandwidth": payload["design"]["bandwidth"], "train": train, "eval": evaluation}, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
