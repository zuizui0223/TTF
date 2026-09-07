#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
from pathlib import Path

from ttf.calibration import CalibrationCell, qualify_calibration

EXPECTED_CELLS = {(0.0, 0.5), (0.0, 1.0), (0.0, 2.0), (0.0, 3.0), (1.0, 2.0)}
N_LAYOUTS = 12


def main() -> int:
    parser = argparse.ArgumentParser(description="Aggregate Queensland topology qualification across all feasible table layouts.")
    parser.add_argument("--input-dir", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    grouped: dict[int, list[CalibrationCell]] = {i: [] for i in range(N_LAYOUTS)}
    seen: dict[int, set[tuple[float, float]]] = {i: set() for i in range(N_LAYOUTS)}
    designs: dict[int, dict] = {}
    for path in sorted(args.input_dir.rglob("*.json")):
        payload = json.loads(path.read_text())
        if payload.get("schema") != "ttf_queensland33_topology_cell_v0.1":
            continue
        layout = int(payload["layout_index"])
        if layout not in grouped:
            raise RuntimeError(f"unexpected layout index {layout}")
        cell = CalibrationCell(**payload["cell"])
        key = (float(cell.shared_fraction), float(cell.amplitude))
        if key in seen[layout]:
            raise RuntimeError(f"duplicate cell for layout {layout}: {key}")
        seen[layout].add(key)
        grouped[layout].append(cell)
        designs.setdefault(layout, payload["design"])

    per_layout = {}
    max_type1 = -1.0
    min_power = 2.0
    for layout in range(N_LAYOUTS):
        if seen[layout] != EXPECTED_CELLS:
            raise RuntimeError(f"layout {layout} grid drift: {sorted(seen[layout])}")
        cells = sorted(grouped[layout], key=lambda c: (c.shared_fraction, c.amplitude))
        report = qualify_calibration(cells, moderate_amplitude=2.0, type1_ceiling=0.10, power_floor=0.80)
        max_type1 = max(max_type1, float(report.max_zero_shared_rejection))
        min_power = min(min_power, float(report.full_shared_moderate_power))
        per_layout[str(layout)] = {
            "cells": [cell.to_dict() for cell in cells],
            "point_qualification": report.to_dict(),
            "design": designs[layout],
        }

    passed = bool(max_type1 <= 0.10 and min_power >= 0.80)
    payload = {
        "schema": "ttf_queensland33_topology_qualification_v0.1",
        "qualification": {
            "passed": passed,
            "worst_layout_max_zero_shared_rejection": max_type1,
            "worst_layout_full_shared_A2_power": min_power,
            "type1_ceiling": 0.10,
            "power_floor": 0.80,
            "n_layouts": N_LAYOUTS,
            "worlds_per_cell_per_layout": 50,
            "pilot_only": True,
        },
        "selection_rule": {
            "all_layouts_must_pass": True,
            "if_pass": "validate on an independent layout-objective seed and higher world count before empirical genetic outcomes",
            "if_fail": "do not open empirical genetic outcomes; record whether failure is type-I, power, or both",
            "genetic_outcomes_used": False,
            "named_Mary_Brisbane_boundary_used": False,
        },
        "layouts": per_layout,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    print(json.dumps(payload["qualification"], sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
