#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

from ttf.calibration import CalibrationCell, qualify_calibration
from ttf.precision import qualify_calibration_precision

MANDATORY_CELLS = {
    (0.0, 0.5),
    (0.0, 1.0),
    (0.0, 2.0),
    (0.0, 3.0),
    (1.0, 2.0),
}


def canonical(value: object) -> str:
    return json.dumps(value, sort_keys=True, separators=(",", ":"))


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Aggregate the five mandatory high-precision TTF Gate-I geometry cells."
    )
    parser.add_argument("--input-dir", type=Path, required=True)
    parser.add_argument("--moderate-amplitude", type=float, default=2.0)
    parser.add_argument("--type1-upper-ceiling", type=float, default=0.10)
    parser.add_argument("--power-lower-floor", type=float, default=0.80)
    parser.add_argument("--gate", default="I-A external non-flower geometry stress")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    paths = sorted(args.input_dir.glob("cell-*.json"))
    if not paths:
        raise RuntimeError("no Gate-I geometry cells found")

    cells: list[CalibrationCell] = []
    payloads: list[dict] = []
    seen: set[tuple[float, float]] = set()
    for path in paths:
        payload = json.loads(path.read_text())
        if payload.get("schema") != "ttf_gate_i_geometry_cell_v0.1":
            raise RuntimeError(f"unexpected Gate-I cell schema in {path}")
        cell = CalibrationCell(**payload["cell"])
        key = (float(cell.shared_fraction), float(cell.amplitude))
        if key in seen:
            raise RuntimeError(f"duplicate Gate-I cell {key}")
        seen.add(key)
        cells.append(cell)
        payloads.append(payload)

    if seen != MANDATORY_CELLS:
        missing = sorted(MANDATORY_CELLS - seen)
        extra = sorted(seen - MANDATORY_CELLS)
        raise RuntimeError(f"Gate-I precision cell set drifted; missing={missing}, extra={extra}")

    reference = payloads[0]
    for section in ("geometry", "split"):
        expected = canonical(reference[section])
        for payload in payloads[1:]:
            if canonical(payload[section]) != expected:
                raise RuntimeError(f"Gate-I {section} drift across precision cells")

    invariant_keys = (
        "replicates",
        "resamples",
        "eval_fraction",
        "k",
        "max_distance",
        "bandwidth",
        "prior_strength",
        "prior_mean",
        "segment_points",
        "noise_sd",
        "transition_width",
        "alpha",
        "seed",
        "inference",
    )
    common: dict[str, object] = {}
    for key in invariant_keys:
        values = {canonical(payload["config"][key]) for payload in payloads}
        if len(values) != 1:
            raise RuntimeError(f"Gate-I config drift for {key}: {sorted(values)}")
        common[key] = reference["config"][key]

    if int(common["replicates"]) != 500:
        raise RuntimeError("claim-bearing Gate-I precision requires exactly 500 worlds per cell")
    if int(common["resamples"]) != 1999:
        raise RuntimeError("claim-bearing Gate-I precision requires exactly 1,999 bootstrap resamples")
    if not np.isclose(float(common["alpha"]), 0.05):
        raise RuntimeError("claim-bearing Gate-I precision requires alpha=0.05")
    if not np.isclose(float(args.moderate_amplitude), 2.0):
        raise RuntimeError("the prospectively frozen moderate amplitude is 2.0")

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
    cells.sort(key=lambda cell: (cell.shared_fraction, cell.amplitude))

    payload = {
        "schema": "ttf_gate_i_geometry_precision_v0.1",
        "gate": args.gate,
        "qualification_scope": "external_non_flower_geometry_high_precision",
        "geometry": reference["geometry"],
        "split": reference["split"],
        "design": {
            **common,
            "mandatory_cells": [list(cell) for cell in sorted(MANDATORY_CELLS)],
            "moderate_amplitude": args.moderate_amplitude,
            "type1_upper_ceiling": args.type1_upper_ceiling,
            "power_lower_floor": args.power_lower_floor,
            "confidence_interval": "two-sided Wilson score 95%",
        },
        "cells": [cell.to_dict() for cell in cells],
        "point_qualification": point.to_dict(),
        "precision_qualification": precision.to_dict(),
        "gate_i_a_pass": bool(precision.passed),
        "gate_i_b_required": True,
        "overall_gate_i_pass": False,
        "claim_ready": False,
        "claim_boundary": (
            "This run can qualify Gate I-A, the independent non-flower geometry stress. "
            "It cannot qualify Gate I overall: Gate I-B must separately pass on the "
            "prospectively frozen intended empirical TTF sampling geometry."
        ),
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    print(json.dumps(payload["precision_qualification"], sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
