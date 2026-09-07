#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
from pathlib import Path

from ttf.calibration import CalibrationCell
from ttf.precision import qualify_calibration_precision


def main() -> int:
    parser = argparse.ArgumentParser(description="Aggregate high-precision TTF qualification cells.")
    parser.add_argument("--input-dir", type=Path, required=True)
    parser.add_argument("--moderate-amplitude", type=float, default=2.0)
    parser.add_argument("--type1-upper-ceiling", type=float, default=0.10)
    parser.add_argument("--power-lower-floor", type=float, default=0.80)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    paths = sorted(args.input_dir.glob("cell-*.json"))
    if not paths:
        raise RuntimeError("no precision calibration cells found")

    cells: list[CalibrationCell] = []
    configs: list[dict] = []
    seen: set[tuple[float, float]] = set()
    for path in paths:
        payload = json.loads(path.read_text())
        if payload.get("schema") != "ttf_calibration_cell_v0.1":
            raise RuntimeError(f"unexpected cell schema in {path}")
        cell = CalibrationCell(**payload["cell"])
        key = (float(cell.shared_fraction), float(cell.amplitude))
        if key in seen:
            raise RuntimeError(f"duplicate cell {key}")
        seen.add(key)
        cells.append(cell)
        configs.append(payload["config"])

    invariant_keys = [
        "replicates",
        "permutations",
        "species",
        "records_per_species",
        "bandwidth",
        "alpha",
        "seed",
        "inference",
    ]
    common = {}
    for key in invariant_keys:
        values = {json.dumps(config[key], sort_keys=True) for config in configs}
        if len(values) != 1:
            raise RuntimeError(f"config drift for {key}: {sorted(values)}")
        common[key] = configs[0][key]

    report = qualify_calibration_precision(
        cells,
        moderate_amplitude=args.moderate_amplitude,
        type1_upper_ceiling=args.type1_upper_ceiling,
        power_lower_floor=args.power_lower_floor,
    )
    cells.sort(key=lambda cell: (cell.shared_fraction, cell.amplitude))
    payload = {
        "schema": "ttf_precision_qualification_v0.2",
        "design": {
            **common,
            "moderate_amplitude": args.moderate_amplitude,
            "type1_upper_ceiling": args.type1_upper_ceiling,
            "power_lower_floor": args.power_lower_floor,
            "confidence_interval": "two-sided Wilson score 95%",
        },
        "cells": [cell.to_dict() for cell in cells],
        "qualification": report.to_dict(),
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    print(json.dumps(payload["qualification"], sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
