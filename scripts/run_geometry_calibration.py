#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import numpy as np

from ttf.calibration import (
    qualify_calibration,
    run_geometry_calibration,
    seed_for,
)
from ttf.core import split_species
from ttf.geometry import SpeciesGeometry, geometry_fingerprint
from ttf.precision import qualify_calibration_precision


def floats(text: str) -> list[float]:
    values = [float(part.strip()) for part in text.split(",") if part.strip()]
    if not values:
        raise argparse.ArgumentTypeError("expected comma-separated numeric values")
    return values


def columns(text: str) -> list[str]:
    values = [part.strip() for part in text.split(",") if part.strip()]
    if not values:
        raise argparse.ArgumentTypeError("expected comma-separated column names")
    return values


def load_geometry(
    path: Path,
    *,
    species_column: str,
    coordinate_columns: list[str],
    block_column: str | None,
    min_records: int,
) -> tuple[tuple[SpeciesGeometry, ...], tuple[str, ...]]:
    if min_records < 2:
        raise ValueError("min_records must be >= 2")
    coordinates: dict[str, list[list[float]]] = {}
    blocks: dict[str, list[str]] = {}
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        fieldnames = set(reader.fieldnames or ())
        required = {species_column, *coordinate_columns}
        if block_column is not None:
            required.add(block_column)
        missing = required - fieldnames
        if missing:
            raise ValueError(f"missing CSV columns: {sorted(missing)}")
        for row_number, row in enumerate(reader, start=2):
            species = str(row[species_column]).strip()
            if not species:
                raise ValueError(f"empty species at CSV row {row_number}")
            try:
                point = [float(row[name]) for name in coordinate_columns]
            except (TypeError, ValueError) as exc:
                raise ValueError(
                    f"non-numeric coordinate at CSV row {row_number}"
                ) from exc
            if not np.isfinite(point).all():
                raise ValueError(f"non-finite coordinate at CSV row {row_number}")
            coordinates.setdefault(species, []).append(point)
            if block_column is not None:
                blocks.setdefault(species, []).append(str(row[block_column]))

    excluded = tuple(
        sorted(species for species, rows in coordinates.items() if len(rows) < min_records)
    )
    geometries: list[SpeciesGeometry] = []
    for species in sorted(coordinates):
        rows = coordinates[species]
        if len(rows) < min_records:
            continue
        geometries.append(
            SpeciesGeometry(
                species=species,
                coordinates=np.asarray(rows, dtype=float),
                blocks=(
                    None
                    if block_column is None
                    else np.asarray(blocks[species], dtype=object)
                ),
            )
        )
    return tuple(geometries), excluded


def main() -> int:
    parser = argparse.ArgumentParser(
        description=(
            "Run TTF Gate-I semi-synthetic qualification on fixed empirical "
            "sampling geometry. CSV trait columns, if present, are ignored."
        )
    )
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--species-column", default="species")
    parser.add_argument("--coordinate-columns", type=columns, default=["x", "y"])
    parser.add_argument("--block-column", default=None)
    parser.add_argument("--min-records", type=int, default=8)
    parser.add_argument("--shared-fractions", type=floats, default=[0.0, 1.0])
    parser.add_argument("--amplitudes", type=floats, default=[0.5, 1.0, 2.0, 3.0])
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
    parser.add_argument("--moderate-amplitude", type=float, default=2.0)
    parser.add_argument("--type1-upper-ceiling", type=float, default=0.10)
    parser.add_argument("--power-lower-floor", type=float, default=0.80)
    parser.add_argument("--alpha", type=float, default=0.05)
    parser.add_argument("--seed", type=int, default=20260907)
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
        raise ValueError(
            "fewer than six held-out species after filtering; Gate-I inference is undefined"
        )

    cells = run_geometry_calibration(
        geometries,
        shared_fractions=args.shared_fractions,
        amplitudes=args.amplitudes,
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
    )
    point = qualify_calibration(
        cells,
        moderate_amplitude=args.moderate_amplitude,
        type1_ceiling=args.type1_upper_ceiling,
        power_floor=args.power_lower_floor,
    )
    precision = qualify_calibration_precision(
        cells,
        moderate_amplitude=args.moderate_amplitude,
        type1_upper_ceiling=args.type1_upper_ceiling,
        power_lower_floor=args.power_lower_floor,
    )

    record_counts = {item.species: len(item.coordinates) for item in geometries}
    payload = {
        "schema": "ttf_gate_i_geometry_calibration_v0.1",
        "geometry": {
            "input": str(args.input),
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
            "shared_fractions": args.shared_fractions,
            "amplitudes": args.amplitudes,
            "replicates": args.replicates,
            "resamples": args.resamples,
            "k": args.k,
            "max_distance": args.max_distance,
            "bandwidth": args.bandwidth,
            "prior_strength": args.prior_strength,
            "prior_mean": args.prior_mean,
            "segment_points": args.segment_points,
            "noise_sd": args.noise_sd,
            "transition_width": args.transition_width,
            "moderate_amplitude": args.moderate_amplitude,
            "type1_upper_ceiling": args.type1_upper_ceiling,
            "power_lower_floor": args.power_lower_floor,
            "alpha": args.alpha,
            "seed": args.seed,
            "inference": "heldout_species_bootstrap",
        },
        "cells": [cell.to_dict() for cell in cells],
        "point_qualification": point.to_dict(),
        "precision_qualification": precision.to_dict(),
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    print(json.dumps(payload["precision_qualification"], sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
