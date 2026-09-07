#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

from ttf.calibration import CalibrationCell, qualify_calibration


def wilson_interval(successes: int, trials: int, z: float = 1.959963984540054) -> tuple[float, float]:
    if trials < 1 or not 0 <= successes <= trials:
        raise ValueError("invalid binomial counts")
    p = successes / trials
    z2 = z * z
    denom = 1.0 + z2 / trials
    center = (p + z2 / (2.0 * trials)) / denom
    half = z * math.sqrt((p * (1.0 - p) + z2 / (4.0 * trials)) / trials) / denom
    return max(0.0, center - half), min(1.0, center + half)


def main() -> int:
    parser = argparse.ArgumentParser(description="Aggregate TTF calibration-cell artifacts.")
    parser.add_argument("--input-dir", type=Path, required=True)
    parser.add_argument("--moderate-amplitude", type=float, default=2.0)
    parser.add_argument("--type1-ceiling", type=float, default=0.10)
    parser.add_argument("--power-floor", type=float, default=0.80)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    paths = sorted(args.input_dir.glob("cell-*.json"))
    if not paths:
        raise RuntimeError("no calibration cell files found")

    cells: list[CalibrationCell] = []
    uncertainty: list[dict[str, float | int]] = []
    seen: set[tuple[float, float]] = set()
    configs: list[dict] = []
    for path in paths:
        payload = json.loads(path.read_text())
        if payload.get("schema") != "ttf_calibration_cell_v0.1":
            raise RuntimeError(f"unexpected schema in {path}")
        cell = CalibrationCell(**payload["cell"])
        key = (float(cell.shared_fraction), float(cell.amplitude))
        if key in seen:
            raise RuntimeError(f"duplicate calibration cell {key}")
        seen.add(key)
        cells.append(cell)
        configs.append(payload["config"])
        rejected = int(round(cell.rejection_rate * cell.n_replicates))
        lo, hi = wilson_interval(rejected, cell.n_replicates)
        uncertainty.append(
            {
                "shared_fraction": cell.shared_fraction,
                "amplitude": cell.amplitude,
                "rejections": rejected,
                "replicates": cell.n_replicates,
                "rejection_rate": cell.rejection_rate,
                "wilson95_low": lo,
                "wilson95_high": hi,
            }
        )

    report = qualify_calibration(
        cells,
        moderate_amplitude=args.moderate_amplitude,
        type1_ceiling=args.type1_ceiling,
        power_floor=args.power_floor,
    )
    cells.sort(key=lambda c: (c.shared_fraction, c.amplitude))
    uncertainty.sort(key=lambda row: (float(row["shared_fraction"]), float(row["amplitude"])))

    invariant_keys = [
        "replicates",
        "permutations",
        "species",
        "records_per_species",
        "bandwidth",
        "alpha",
        "seed",
    ]
    common = {}
    for key in invariant_keys:
        values = {json.dumps(cfg[key], sort_keys=True) for cfg in configs}
        if len(values) != 1:
            raise RuntimeError(f"calibration config drift for {key}: {sorted(values)}")
        common[key] = configs[0][key]

    payload = {
        "schema": "ttf_qualification_pilot_v0.1",
        "design": {
            **common,
            "shared_fractions": sorted({c.shared_fraction for c in cells}),
            "amplitudes": sorted({c.amplitude for c in cells}),
            "moderate_amplitude": args.moderate_amplitude,
            "type1_ceiling": args.type1_ceiling,
            "power_floor": args.power_floor,
            "uncertainty": "Wilson score 95% interval; provisional pilot only",
        },
        "cells": [cell.to_dict() for cell in cells],
        "binomial_uncertainty": uncertainty,
        "qualification": report.to_dict(),
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    print(json.dumps(payload["qualification"], sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
