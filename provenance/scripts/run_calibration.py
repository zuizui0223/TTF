#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
from pathlib import Path

from ttf.calibration import qualify_calibration, run_calibration


def floats(text: str) -> list[float]:
    values = [float(part.strip()) for part in text.split(",") if part.strip()]
    if not values:
        raise argparse.ArgumentTypeError("expected comma-separated numeric values")
    return values


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Run the prospective TTF sharedness x amplitude calibration grid."
    )
    parser.add_argument("--shared-fractions", type=floats, default=[0.0, 0.25, 0.5, 1.0])
    parser.add_argument("--amplitudes", type=floats, default=[0.0, 1.0, 2.0, 3.0])
    parser.add_argument("--replicates", type=int, default=100)
    parser.add_argument("--resamples", "--permutations", dest="resamples", type=int, default=1999)
    parser.add_argument("--species", type=int, default=40)
    parser.add_argument("--records-per-species", type=int, default=60)
    parser.add_argument("--bandwidth", type=float, default=0.2)
    parser.add_argument("--moderate-amplitude", type=float, default=2.0)
    parser.add_argument("--alpha", type=float, default=0.05)
    parser.add_argument("--seed", type=int, default=20260907)
    parser.add_argument(
        "--inference",
        choices=["heldout_species_bootstrap", "trait_permutation"],
        default="heldout_species_bootstrap",
        help="Primary sharedness inference defaults to held-out species bootstrap; trait permutation is diagnostic.",
    )
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    cells = run_calibration(
        shared_fractions=args.shared_fractions,
        amplitudes=args.amplitudes,
        n_replicates=args.replicates,
        n_permutations=args.resamples,
        n_species=args.species,
        records_per_species=args.records_per_species,
        bandwidth=args.bandwidth,
        alpha=args.alpha,
        seed=args.seed,
        inference=args.inference,
    )
    qualification = qualify_calibration(
        cells,
        moderate_amplitude=args.moderate_amplitude,
    )
    payload = {
        "schema": "ttf_calibration_v0.2",
        "config": {
            "shared_fractions": args.shared_fractions,
            "amplitudes": args.amplitudes,
            "replicates": args.replicates,
            "resamples": args.resamples,
            "species": args.species,
            "records_per_species": args.records_per_species,
            "bandwidth": args.bandwidth,
            "moderate_amplitude": args.moderate_amplitude,
            "alpha": args.alpha,
            "seed": args.seed,
            "inference": args.inference,
        },
        "cells": [cell.to_dict() for cell in cells],
        "qualification": qualification.to_dict(),
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    print(json.dumps(payload["qualification"], sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
