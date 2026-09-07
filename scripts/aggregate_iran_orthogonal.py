#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
from pathlib import Path

from ttf.calibration import CalibrationCell, qualify_calibration

VARIANTS = ("orthogonal_opportunity", "orthogonal_geometry")
EXPECTED = {(0.0, 0.5), (0.0, 1.0), (0.0, 2.0), (0.0, 3.0), (1.0, 2.0)}


def main() -> int:
    parser = argparse.ArgumentParser(description="Aggregate Iranian predictor-orthogonalized TTF cells.")
    parser.add_argument("--input-dir", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    grouped = {variant: [] for variant in VARIANTS}
    designs = {variant: [] for variant in VARIANTS}
    diagnostics = {variant: [] for variant in VARIANTS}
    seen = {variant: set() for variant in VARIANTS}

    for path in sorted(args.input_dir.rglob("*.json")):
        payload = json.loads(path.read_text())
        if payload.get("schema") != "ttf_iran_orthogonal_cell_v0.1":
            continue
        variant = str(payload["variant"])
        if variant not in grouped:
            raise RuntimeError(f"unexpected orthogonal variant {variant}")
        cell = CalibrationCell(**payload["cell"])
        key = (float(cell.shared_fraction), float(cell.amplitude))
        if key in seen[variant]:
            raise RuntimeError(f"duplicate cell for {variant}: {key}")
        seen[variant].add(key)
        grouped[variant].append(cell)
        designs[variant].append(payload["design"])
        diagnostics[variant].append({
            "shared_fraction": float(cell.shared_fraction),
            "amplitude": float(cell.amplitude),
            **payload["diagnostic_means"],
        })

    output_variants = {}
    eligible = []
    for variant in VARIANTS:
        if seen[variant] != EXPECTED:
            raise RuntimeError(f"{variant} grid drift: {sorted(seen[variant])}")
        if len({json.dumps(item, sort_keys=True) for item in designs[variant]}) != 1:
            raise RuntimeError(f"design drift across {variant} cells")
        cells = sorted(grouped[variant], key=lambda item: (item.shared_fraction, item.amplitude))
        report = qualify_calibration(cells, moderate_amplitude=2.0, type1_ceiling=0.10, power_floor=0.80)
        if report.passed:
            eligible.append(variant)
        output_variants[variant] = {
            "cells": [cell.to_dict() for cell in cells],
            "diagnostics": sorted(diagnostics[variant], key=lambda row: (row["shared_fraction"], row["amplitude"])),
            "point_qualification": report.to_dict(),
            "design": designs[variant][0],
        }

    payload = {
        "schema": "ttf_iran_orthogonal_qualification_v0.1",
        "variants": output_variants,
        "selection_rule": {
            "empirical_genetic_outcomes_used": False,
            "eligible_for_independent_seed_validation": "point type-I <= 0.10 and power >= 0.80 on all frozen cells",
            "eligible_variants": eligible,
            "if_multiple_pass": "validate all passing variants on an independent world/bootstrap seed before choosing",
            "if_none_pass": "do not open empirical SNP outcomes; redesign the estimator/null",
        },
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    print(json.dumps({
        "eligible": eligible,
        **{
            variant: {
                "passed": output_variants[variant]["point_qualification"]["passed"],
                "max_type1": output_variants[variant]["point_qualification"]["max_zero_shared_rejection"],
                "power": output_variants[variant]["point_qualification"]["full_shared_moderate_power"],
            }
            for variant in VARIANTS
        },
    }, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
