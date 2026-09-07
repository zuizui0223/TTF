#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
from pathlib import Path

from ttf.calibration import CalibrationCell, qualify_calibration
from ttf.precision import qualify_calibration_precision

EXPECTED = {(0.0, 0.5), (0.0, 1.0), (0.0, 2.0), (0.0, 3.0), (1.0, 2.0)}


def main() -> int:
    parser = argparse.ArgumentParser(description="Aggregate Iranian actual-geometry TTF calibration cells.")
    parser.add_argument("--input-dir", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    paths = sorted(args.input_dir.rglob("cell-*.json"))
    if not paths:
        raise RuntimeError("no Iranian geometry calibration cells found")

    cells: list[CalibrationCell] = []
    designs: list[dict] = []
    seen: set[tuple[float, float]] = set()
    for path in paths:
        payload = json.loads(path.read_text())
        if payload.get("schema") != "ttf_iran_geometry_calibration_cell_v0.1":
            continue
        cell = CalibrationCell(**payload["cell"])
        key = (float(cell.shared_fraction), float(cell.amplitude))
        if key in seen:
            raise RuntimeError(f"duplicate Iranian geometry cell {key}")
        seen.add(key)
        cells.append(cell)
        designs.append(payload["design"])

    if seen != EXPECTED:
        raise RuntimeError(f"Iranian geometry calibration cells differ from frozen grid: {sorted(seen)}")

    invariant = [
        "inference", "geometry_schema", "coordinate_projection", "bandwidth_rule",
        "bandwidth", "k", "noise_sd", "transition_width", "min_shared_crossing_species",
        "split_seed", "master_seed", "train_species", "eval_species", "bootstrap_resamples",
        "outcome_values_used", "known_boundary_labels_used",
    ]
    common = {}
    for key in invariant:
        serialized = {json.dumps(design[key], sort_keys=True) for design in designs}
        if len(serialized) != 1:
            raise RuntimeError(f"design drift for {key}: {serialized}")
        common[key] = designs[0][key]

    cells.sort(key=lambda item: (item.shared_fraction, item.amplitude))
    point = qualify_calibration(cells, moderate_amplitude=2.0, type1_ceiling=0.10, power_floor=0.80)
    precision = qualify_calibration_precision(
        cells,
        moderate_amplitude=2.0,
        type1_upper_ceiling=0.10,
        power_lower_floor=0.80,
    )
    payload = {
        "schema": "ttf_iran_geometry_qualification_v0.1",
        "design": {
            **common,
            "replicates_per_cell": cells[0].n_replicates,
            "alpha": cells[0].alpha,
            "moderate_amplitude": 2.0,
            "point_type1_ceiling": 0.10,
            "point_power_floor": 0.80,
            "precision_type1_upper_ceiling": 0.10,
            "precision_power_lower_floor": 0.80,
        },
        "cells": [cell.to_dict() for cell in cells],
        "point_qualification": point.to_dict(),
        "precision_qualification": precision.to_dict(),
        "interpretation": {
            "point_pass_is_pilot_only": True,
            "empirical_genetic_outcomes_may_be_opened_only_after_precision_pass": True,
        },
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    print(json.dumps({
        "point_pass": point.passed,
        "precision_pass": precision.passed,
        "max_type1": point.max_zero_shared_rejection,
        "power": point.full_shared_moderate_power,
        "max_type1_upper95": precision.max_zero_shared_upper95,
        "power_lower95": precision.full_shared_moderate_lower95,
    }, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
