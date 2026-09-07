#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
from pathlib import Path

from ttf.calibration import CalibrationCell
from ttf.precision import qualify_calibration_precision

PANELS = ("v03", "v04")
VARIANTS = (
    "prediction_propensity_only",
    "prediction_length_and_propensity",
    "target_propensity_only",
    "two_sided_propensity_only",
    "v04_two_sided_length_and_propensity",
)
PREFERENCE = (
    "prediction_propensity_only",
    "prediction_length_and_propensity",
    "target_propensity_only",
    "two_sided_propensity_only",
    "v04_two_sided_length_and_propensity",
)
MANDATORY = {(0.0,0.5),(0.0,1.0),(0.0,2.0),(0.0,3.0),(1.0,2.0)}


def main() -> int:
    parser = argparse.ArgumentParser(description="Aggregate v0.5 one-sided geometry diagnostics across both failed external panels.")
    parser.add_argument("--input-dir", type=Path, required=True)
    parser.add_argument("--rule", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    rule = json.loads(args.rule.read_text())
    if rule.get("status") != "frozen_before_v05_diagnostic_outcomes":
        raise RuntimeError("v0.5 rule was not prospectively frozen")
    if tuple(rule.get("candidate_preference_order", ())) != PREFERENCE:
        raise RuntimeError("candidate preference order drifted")
    req = rule["precision_requirements"]
    if req["worlds_per_cell"] != 500 or req["bootstrap_resamples"] != 1999:
        raise RuntimeError("precision execution drifted")

    paths = sorted(args.input_dir.glob("cell-*.json"))
    if len(paths) != 10:
        raise RuntimeError(f"expected 10 diagnostic cells, found {len(paths)}")
    data = {panel: {} for panel in PANELS}
    for path in paths:
        payload = json.loads(path.read_text())
        if payload.get("schema") != "ttf_v05_one_sided_diagnostic_cell_v0.1":
            raise RuntimeError(f"unexpected schema in {path}")
        if payload.get("status") != "development_only" or payload.get("baseline_reproduced_exactly") is not True:
            raise RuntimeError(f"development/baseline firewall failed in {path}")
        panel = payload["panel"]
        if panel not in PANELS:
            raise RuntimeError(f"unexpected panel {panel}")
        key=(float(payload["shared_fraction"]), float(payload["amplitude"]))
        if key in data[panel]:
            raise RuntimeError(f"duplicate {panel} cell {key}")
        data[panel][key]=payload
    for panel in PANELS:
        if set(data[panel]) != MANDATORY:
            raise RuntimeError(f"{panel} mandatory cells drifted")

    candidates = {}
    for variant in VARIANTS:
        panel_reports = {}
        all_pass = True
        for panel in PANELS:
            cells=[]
            for key in sorted(MANDATORY):
                cells.append(CalibrationCell(**data[panel][key]["variants"][variant]))
            precision = qualify_calibration_precision(
                cells,
                moderate_amplitude=2.0,
                type1_upper_ceiling=0.10,
                power_lower_floor=0.80,
            )
            panel_reports[panel] = {
                "cells": [cell.to_dict() for cell in cells],
                "precision": precision.to_dict(),
                "passed": bool(precision.passed),
            }
            all_pass = all_pass and bool(precision.passed)
        candidates[variant] = {
            "panels": panel_reports,
            "passes_both_development_panels": all_pass,
        }

    selected = None
    for variant in PREFERENCE:
        if candidates[variant]["passes_both_development_panels"]:
            selected=variant
            break

    out={
        "schema":"ttf_v05_one_sided_diagnostic_v0.1",
        "status":"development_only_after_v03_v04_external_failures",
        "rule":str(args.rule),
        "candidate_preference_order":list(PREFERENCE),
        "candidates":candidates,
        "selected_candidate":selected,
        "selection_rule_applied_mechanically":True,
        "claim_ready":False,
        "failed_panels_cannot_qualify_v05":True,
        "confirmation_requires_predeclared_birds_and_butterflies":True,
        "rgfca_reserve_opened":False,
        "interpretation_boundary":"Both failed external panels are development-only. A selected candidate, if any, must be locked before the predeclared bird and butterfly confirmation panels are opened; both confirmations must pass Wilson precision before the RGFCA reserve can be authorized."
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(out, indent=2, sort_keys=True)+"\n")
    print(json.dumps({
        "selected_candidate":selected,
        "passes":{v:candidates[v]["passes_both_development_panels"] for v in VARIANTS},
        "panel_precision":{
            v:{p:candidates[v]["panels"][p]["precision"] for p in PANELS}
            for v in VARIANTS
        },
    }, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
