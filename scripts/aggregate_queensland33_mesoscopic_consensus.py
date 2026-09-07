#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np


def main() -> int:
    parser = argparse.ArgumentParser(description="Aggregate Queensland mesoscopic-consensus qualification cells.")
    parser.add_argument("--input-dir", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--type1-ceiling", type=float, default=0.10)
    parser.add_argument("--power-floor", type=float, default=0.80)
    args = parser.parse_args()

    payloads = []
    for path in sorted(args.input_dir.glob("*.json")):
        payload = json.loads(path.read_text())
        if payload.get("schema") == "ttf_queensland33_mesoscopic_consensus_cell_v0.3b":
            payloads.append(payload)
    if not payloads:
        raise RuntimeError("no v0.3b consensus calibration cells found")

    cells = {
        (
            int(payload["layout_index"]),
            float(payload["cell"]["shared_fraction"]),
            float(payload["cell"]["amplitude"]),
        ): payload
        for payload in payloads
    }
    layouts = set(range(12))
    if {key[0] for key in cells} != layouts:
        raise RuntimeError("v0.3b requires all 12 frozen layouts")

    private_amplitudes = (0.5, 1.0, 2.0, 3.0, 4.0)
    shared_fractions = (0.25, 0.5, 1.0)
    for layout in layouts:
        for amplitude in private_amplitudes:
            if (layout, 0.0, amplitude) not in cells:
                raise RuntimeError(f"missing private cell layout={layout} A={amplitude}")
        for shared in shared_fractions:
            if (layout, shared, 2.0) not in cells:
                raise RuntimeError(f"missing shared cell layout={layout} pi={shared}")

    per_layout: dict[str, dict] = {}
    for layout in sorted(layouts):
        max_type1 = max(
            float(cells[(layout, 0.0, amplitude)]["cell"]["rejection_rate"])
            for amplitude in private_amplitudes
        )
        power = float(cells[(layout, 1.0, 2.0)]["cell"]["rejection_rate"])
        per_layout[str(layout)] = {
            "max_zero_shared_rejection": max_type1,
            "full_shared_A2_power": power,
            "passed": bool(max_type1 <= float(args.type1_ceiling) and power >= float(args.power_floor)),
        }

    worst_type1 = max(value["max_zero_shared_rejection"] for value in per_layout.values())
    worst_power = min(value["full_shared_A2_power"] for value in per_layout.values())
    passed = all(value["passed"] for value in per_layout.values())

    response = {}
    for shared in (0.0, 0.25, 0.5, 1.0):
        subset = [
            payload["cell"]
            for (layout, sf, amplitude), payload in cells.items()
            if np.isclose(sf, shared) and np.isclose(amplitude, 2.0)
        ]
        response[str(shared)] = {
            "mean_rejection_rate_across_layouts": float(np.mean([x["rejection_rate"] for x in subset])),
            "min_rejection_rate_across_layouts": float(np.min([x["rejection_rate"] for x in subset])),
            "max_rejection_rate_across_layouts": float(np.max([x["rejection_rate"] for x in subset])),
            "mean_statistic_across_layouts": float(np.mean([x["mean_statistic"] for x in subset])),
        }

    if passed:
        status = "PASS_PILOT"
    elif worst_type1 > float(args.type1_ceiling) and worst_power < float(args.power_floor):
        status = "FAIL_TYPE1_AND_POWER"
    elif worst_type1 > float(args.type1_ceiling):
        status = "FAIL_TYPE1"
    else:
        status = "FAIL_POWER"

    output = {
        "schema": "ttf_queensland33_mesoscopic_consensus_qualification_v0.3b",
        "design": {
            "candidate_registered_before_v0.3a_result": True,
            "n_layouts": 12,
            "n_eligible_species": 21,
            "n_candidate_interbasin_cuts": 7,
            "split": "10 train / 11 evaluation",
            "k": 2,
            "worlds_per_cell_per_layout": int(payloads[0]["cell"]["n_replicates"]),
            "bootstrap_resamples": int(payloads[0]["design"]["bootstrap_resamples"]),
            "noise_sd": float(payloads[0]["design"]["noise_sd"]),
            "transition_width": float(payloads[0]["design"]["transition_width"]),
            "prior_strength": float(payloads[0]["design"]["prior_strength"]),
            "evidence_temperature": float(payloads[0]["design"]["evidence_temperature"]),
            "exact_frozen_layout_envelope": True,
            "genetic_outcomes_used": False,
            "named_Mary_Brisbane_boundary_used": False,
        },
        "prospective_gate": {
            "type1_ceiling": float(args.type1_ceiling),
            "full_shared_A2_power_floor": float(args.power_floor),
            "all_layouts_must_pass": True,
        },
        "qualification": {
            "passed": bool(passed),
            "pilot_only": True,
            "status": status,
            "worst_layout_full_shared_A2_power": float(worst_power),
            "worst_layout_max_zero_shared_rejection": float(worst_type1),
        },
        "response_to_shared_fraction_A2": response,
        "per_layout": per_layout,
        "interpretation": {
            "null_contract": "Held-out score is a soft-label cross-entropy gain over the same species' geometry-private uniform cut prior. Under a calibrated private posterior its conditional expectation is non-positive for any frozen training consensus.",
            "selection": "The estimator, prior_strength=1, evidence_temperature=1, fresh seed 20260911, and the unchanged 0.10/0.80 Gate were registered before any v0.3a Gate result was inspected.",
            "if_pass": "Confirm on another independent seed and higher world count before opening empirical genetic outcomes.",
            "if_fail": "Retain the result and diagnose posterior calibration or consensus learning without using empirical outcomes.",
        },
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(output, indent=2, sort_keys=True) + "\n")
    print(json.dumps(output["qualification"], sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
