#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

from ttf.precision import wilson_interval


PRIVATE_AMPLITUDES = (0.5, 1.0, 2.0, 3.0, 4.0)


def _interval(cell: dict) -> dict:
    n = int(cell["n_replicates"])
    raw = float(cell["rejection_rate"]) * n
    k = int(round(raw))
    if not np.isclose(raw, k, atol=1e-8, rtol=0.0):
        raise RuntimeError("rejection rate is incompatible with replicate count")
    return wilson_interval(k, n).to_dict()


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Aggregate high-replicate Queensland v0.3c randomization qualification."
    )
    parser.add_argument("--input-dir", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--type1-upper-ceiling", type=float, default=0.10)
    parser.add_argument("--power-lower-floor", type=float, default=0.80)
    args = parser.parse_args()

    payloads = []
    for path in sorted(args.input_dir.glob("*.json")):
        payload = json.loads(path.read_text())
        if payload.get("schema") != "ttf_queensland33_mesoscopic_randomization_cell_v0.3c":
            continue
        payloads.append(payload)
    if not payloads:
        raise RuntimeError("no v0.3c precision cells found")

    cells: dict[tuple[int, float, float], dict] = {}
    for payload in payloads:
        cell = payload["cell"]
        key = (
            int(payload["layout_index"]),
            float(cell["shared_fraction"]),
            float(cell["amplitude"]),
        )
        if key in cells:
            raise RuntimeError(f"duplicate cell {key}")
        cells[key] = payload

    expected = {
        (layout, 0.0, amplitude)
        for layout in range(12)
        for amplitude in PRIVATE_AMPLITUDES
    } | {(layout, 1.0, 2.0) for layout in range(12)}
    if set(cells) != expected:
        missing = sorted(expected - set(cells))
        extra = sorted(set(cells) - expected)
        raise RuntimeError(f"precision cell coverage mismatch; missing={missing}, extra={extra}")

    invariant_keys = (
        "split_seed",
        "k",
        "noise_sd",
        "transition_width",
        "min_shared_crossing_species",
        "prior_strength",
        "alignment_randomizations",
        "master_seed",
    )
    for key in invariant_keys:
        values = {json.dumps(p["design"][key], sort_keys=True) for p in payloads}
        if len(values) != 1:
            raise RuntimeError(f"design drift for {key}: {sorted(values)}")

    replicates = {int(p["cell"]["n_replicates"]) for p in payloads}
    if len(replicates) != 1:
        raise RuntimeError(f"replicate-count drift: {sorted(replicates)}")
    n_worlds = next(iter(replicates))

    per_layout = {}
    all_private_highs = []
    all_power_lows = []
    for layout in range(12):
        private_rows = []
        for amplitude in PRIVATE_AMPLITUDES:
            cell = cells[(layout, 0.0, amplitude)]["cell"]
            interval = _interval(cell)
            all_private_highs.append(float(interval["high"]))
            private_rows.append({
                "amplitude": float(amplitude),
                "rejection_rate": float(cell["rejection_rate"]),
                "wilson95": interval,
                "passed": bool(float(interval["high"]) <= float(args.type1_upper_ceiling)),
            })

        full_cell = cells[(layout, 1.0, 2.0)]["cell"]
        full_interval = _interval(full_cell)
        all_power_lows.append(float(full_interval["low"]))
        layout_type1 = all(row["passed"] for row in private_rows)
        layout_power = float(full_interval["low"]) >= float(args.power_lower_floor)
        per_layout[str(layout)] = {
            "private": private_rows,
            "full_shared_A2": {
                "power": float(full_cell["rejection_rate"]),
                "wilson95": full_interval,
                "passed": bool(layout_power),
            },
            "type1_pass": bool(layout_type1),
            "power_pass": bool(layout_power),
            "passed": bool(layout_type1 and layout_power),
        }

    type1_pass = all(v["type1_pass"] for v in per_layout.values())
    power_pass = all(v["power_pass"] for v in per_layout.values())
    passed = bool(type1_pass and power_pass)
    output = {
        "schema": "ttf_queensland33_mesoscopic_randomization_precision_v0.3c",
        "status": "PASS_PRECISION" if passed else (
            "FAIL_TYPE1_AND_POWER" if not type1_pass and not power_pass else
            "FAIL_TYPE1" if not type1_pass else "FAIL_POWER"
        ),
        "design": {
            "n_layouts": 12,
            "n_eligible_species": 21,
            "n_candidate_interbasin_cuts": 7,
            "split": "10 train / 11 evaluation",
            "worlds_per_cell_per_layout": int(n_worlds),
            "alignment_randomizations": int(payloads[0]["design"]["alignment_randomizations"]),
            "master_seed": int(payloads[0]["design"]["master_seed"]),
            "private_amplitudes": list(PRIVATE_AMPLITUDES),
            "full_shared_positive_control": {"shared_fraction": 1.0, "amplitude": 2.0},
            "confidence_interval": "two-sided Wilson score 95%",
            "method_changed_after_v0.3c_pilot": False,
            "exact_frozen_layout_envelope": True,
            "genetic_outcomes_used": False,
            "named_Mary_Brisbane_boundary_used": False,
        },
        "prospective_precision_gate": {
            "every_private_cell_wilson95_upper_at_or_below": float(args.type1_upper_ceiling),
            "every_layout_full_shared_A2_wilson95_lower_at_or_above": float(args.power_lower_floor),
            "all_layouts_must_pass": True,
        },
        "qualification": {
            "passed": passed,
            "type1_pass": bool(type1_pass),
            "power_pass": bool(power_pass),
            "worst_private_wilson95_upper": float(max(all_private_highs)),
            "worst_full_shared_A2_wilson95_lower": float(min(all_power_lows)),
            "max_private_point_rejection": float(max(
                cells[(layout, 0.0, amplitude)]["cell"]["rejection_rate"]
                for layout in range(12)
                for amplitude in PRIVATE_AMPLITUDES
            )),
            "min_full_shared_A2_point_power": float(min(
                cells[(layout, 1.0, 2.0)]["cell"]["rejection_rate"]
                for layout in range(12)
            )),
        },
        "per_layout": per_layout,
        "interpretation": {
            "if_pass": "The unchanged v0.3c estimator and alignment-randomization null satisfy the predeclared dataset-specific precision Gate on the complete frozen Queensland sampling-layout envelope.",
            "if_fail": "Do not tune the method or open empirical genetic outcomes. Treat the failed bound as a qualification failure and diagnose identifiability using synthetic worlds only."
        },
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(output, indent=2, sort_keys=True) + "\n")
    print(json.dumps(output["qualification"] | {"status": output["status"]}, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
