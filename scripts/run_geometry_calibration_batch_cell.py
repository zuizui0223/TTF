#!/usr/bin/env python3
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

from run_geometry_calibration import columns, load_geometry
from ttf.calibration import seed_for
from ttf.core import split_species
from ttf.geometry import geometry_fingerprint
from ttf.geometry_batch import run_geometry_calibration_batched


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Run one exact batched TTF Gate-I empirical-geometry calibration cell."
    )
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--species-column", default="species")
    parser.add_argument("--coordinate-columns", type=columns, default=["x", "y"])
    parser.add_argument("--block-column", default=None)
    parser.add_argument("--min-records", type=int, default=8)
    parser.add_argument("--shared-fraction", type=float, required=True)
    parser.add_argument("--amplitude", type=float, required=True)
    parser.add_argument("--replicates", type=int, default=500)
    parser.add_argument("--resamples", type=int, default=1999)
    parser.add_argument("--eval-fraction", type=float, default=0.5)
    parser.add_argument("--split-seed", type=int, default=None)
    parser.add_argument("--k", type=int, default=4)
    parser.add_argument("--max-distance", type=float, default=None)
    parser.add_argument("--bandwidth", type=float, required=True)
    parser.add_argument("--prior-strength", type=float, default=0.25)
    parser.add_argument("--prior-mean", type=float, default=0.5)
    parser.add_argument("--segment-points", type=int, default=5)
    parser.add_argument("--noise-sd", type=float, default=0.8)
    parser.add_argument("--transition-width", type=float, default=0.2)
    parser.add_argument("--alpha", type=float, default=0.05)
    parser.add_argument("--seed", type=int, default=20260907)
    parser.add_argument("--world-batch-size", type=int, default=25)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    geometries, excluded = load_geometry(
        args.input,
        species_column=args.species_column,
        coordinate_columns=args.coordinate_columns,
        block_column=args.block_column,
        min_records=args.min_records,
    )
    labels = [item.species for item in geometries]
    use_split_seed = (
        seed_for(args.seed, "fixed_geometry_split")
        if args.split_seed is None
        else int(args.split_seed)
    )
    train, evaluation = split_species(
        labels,
        eval_fraction=args.eval_fraction,
        seed=use_split_seed,
    )
    if len(evaluation) < 6:
        raise ValueError("Gate-I requires at least six held-out species")

    cells = run_geometry_calibration_batched(
        geometries,
        shared_fractions=[args.shared_fraction],
        amplitudes=[args.amplitude],
        n_replicates=args.replicates,
        n_bootstrap=args.resamples,
        train_species=train,
        eval_species=evaluation,
        k=args.k,
        max_distance=args.max_distance,
        bandwidth=args.bandwidth,
        prior_strength=args.prior_strength,
        prior_mean=args.prior_mean,
        segment_points=args.segment_points,
        noise_sd=args.noise_sd,
        transition_width=args.transition_width,
        alpha=args.alpha,
        seed=args.seed,
        world_batch_size=args.world_batch_size,
    )
    if len(cells) != 1:
        raise RuntimeError("one-cell Gate-I runner produced an unexpected grid")

    record_counts = {item.species: len(item.coordinates) for item in geometries}
    payload = {
        "schema": "ttf_gate_i_geometry_batch_cell_v0.1",
        "geometry": {
            "input": str(args.input),
            "input_sha256": hashlib.sha256(args.input.read_bytes()).hexdigest(),
            "fingerprint_sha256": geometry_fingerprint(geometries),
            "coordinate_columns": args.coordinate_columns,
            "species_column": args.species_column,
            "block_column": args.block_column,
            "min_records": args.min_records,
            "species_count": len(geometries),
            "record_count": int(sum(record_counts.values())),
            "records_per_species": record_counts,
            "excluded_below_min_records": list(excluded),
        },
        "split": {
            "split_seed": use_split_seed,
            "train_species": list(train),
            "eval_species": list(evaluation),
        },
        "config": {
            "shared_fraction": float(args.shared_fraction),
            "amplitude": float(args.amplitude),
            "replicates": args.replicates,
            "resamples": args.resamples,
            "eval_fraction": args.eval_fraction,
            "k": args.k,
            "max_distance": args.max_distance,
            "bandwidth": args.bandwidth,
            "prior_strength": args.prior_strength,
            "prior_mean": args.prior_mean,
            "segment_points": args.segment_points,
            "noise_sd": args.noise_sd,
            "transition_width": args.transition_width,
            "alpha": args.alpha,
            "seed": args.seed,
            "inference": "heldout_species_bootstrap",
            "execution": "exact_dense_projection_batched_worlds",
            "world_batch_size": args.world_batch_size,
        },
        "cell": cells[0].to_dict(),
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    print(json.dumps(payload["cell"], sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
