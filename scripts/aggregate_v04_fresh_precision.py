#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

from ttf.calibration import CalibrationCell, qualify_calibration
from ttf.geometry_batch_v04 import (
    V04_LOCKED_CANDIDATE,
    V04_PROPENSITY_DIRECTIONS,
    V04_PROPENSITY_SEED,
    V04_PROPENSITY_TRANSITION_WIDTH,
)
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
    parser = argparse.ArgumentParser(description="Aggregate locked v0.4 fresh external qualification")
    parser.add_argument("--input-dir", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    paths = sorted(args.input_dir.glob("cell-*.json"))
    if len(paths) != 5:
        raise RuntimeError(f"expected exactly five v0.4 cells, got {len(paths)}")
    payloads = [json.loads(path.read_text()) for path in paths]
    for path, payload in zip(paths, payloads):
        if payload.get("schema") != "ttf_v04_fresh_geometry_cell_v0.1":
            raise RuntimeError(f"unexpected v0.4 cell schema in {path}")
        if payload.get("candidate") != V04_LOCKED_CANDIDATE:
            raise RuntimeError("candidate drift across v0.4 cells")
        if payload.get("qualification_role") != "second_fresh_external_confirmatory_geometry":
            raise RuntimeError("v0.4 cell is not second fresh confirmatory geometry")

    reference = payloads[0]
    for section in ("geometry", "split"):
        expected = canonical(reference[section])
        if any(canonical(payload[section]) != expected for payload in payloads[1:]):
            raise RuntimeError(f"v0.4 {section} drift across cells")

    invariant_keys = (
        "replicates", "resamples", "eval_fraction", "k", "bandwidth",
        "prior_strength", "prior_mean", "segment_points", "noise_sd",
        "transition_width", "alpha", "seed", "world_batch_size", "inference",
        "execution", "training_response", "evaluation_score",
        "propensity_directions", "propensity_seed", "propensity_transition_width",
    )
    common: dict[str, object] = {}
    for key in invariant_keys:
        values = {canonical(payload["config"][key]) for payload in payloads}
        if len(values) != 1:
            raise RuntimeError(f"v0.4 config drift for {key}: {sorted(values)}")
        common[key] = reference["config"][key]

    if int(common["replicates"]) != 500 or int(common["resamples"]) != 1999:
        raise RuntimeError("v0.4 precision requires 500 worlds and 1,999 bootstraps per cell")
    if not np.isclose(float(common["alpha"]), 0.05):
        raise RuntimeError("v0.4 qualification requires alpha=0.05")
    if not np.isclose(float(common["prior_mean"]), 0.0):
        raise RuntimeError("v0.4 inherits centered v0.3 training field")
    if common["training_response"] != "v03_length_rank_orthogonalized_turnover":
        raise RuntimeError("v0.4 training response drifted")
    if common["evaluation_score"] != "partial_spearman_control_edge_length_and_random_private_geometry_propensity":
        raise RuntimeError("v0.4 evaluation score drifted")
    if int(common["propensity_directions"]) != V04_PROPENSITY_DIRECTIONS:
        raise RuntimeError("v0.4 propensity direction-count drift")
    if int(common["propensity_seed"]) != V04_PROPENSITY_SEED:
        raise RuntimeError("v0.4 propensity seed drift")
    if not np.isclose(float(common["propensity_transition_width"]), V04_PROPENSITY_TRANSITION_WIDTH):
        raise RuntimeError("v0.4 propensity width drift")

    cells: list[CalibrationCell] = []
    seen = set()
    for payload in payloads:
        cell = CalibrationCell(**payload["cell"])
        key = (float(cell.shared_fraction), float(cell.amplitude))
        if key in seen:
            raise RuntimeError(f"duplicate v0.4 cell {key}")
        seen.add(key)
        cells.append(cell)
    if seen != MANDATORY_CELLS:
        raise RuntimeError(f"v0.4 mandatory cell set drifted: {sorted(seen)}")

    point = qualify_calibration(
        cells,
        moderate_amplitude=2.0,
        type1_ceiling=0.10,
        power_floor=0.80,
    )
    precision = qualify_calibration_precision(
        cells,
        moderate_amplitude=2.0,
        type1_upper_ceiling=0.10,
        power_lower_floor=0.80,
    )
    cells.sort(key=lambda cell: (cell.shared_fraction, cell.amplitude))

    result = {
        "schema": "ttf_v04_fresh_external_precision_v0.1",
        "candidate": V04_LOCKED_CANDIDATE,
        "qualification_scope": "second_fresh_external_animal_geometry_only",
        "fresh_geometry": reference["geometry"],
        "split": reference["split"],
        "design": {
            **common,
            "mandatory_cells": [list(cell) for cell in sorted(MANDATORY_CELLS)],
            "moderate_amplitude": 2.0,
            "type1_upper_ceiling": 0.10,
            "power_lower_floor": 0.80,
            "confidence_interval": "two-sided Wilson score 95%",
        },
        "cells": [cell.to_dict() for cell in cells],
        "point_qualification": point.to_dict(),
        "precision_qualification": precision.to_dict(),
        "v04_fresh_external_pass": bool(precision.passed),
        "claim_ready": False,
        "reserve_gate_permitted": bool(precision.passed),
        "claim_boundary": (
            "A pass qualifies the locked v0.4 evaluation-only geometry partial on this second "
            "prospectively frozen external animal geometry. It does not retroactively alter v0.2/v0.3 "
            "failures and does not itself qualify the intended RGFCA domain."
        ),
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps(result["precision_qualification"], sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
