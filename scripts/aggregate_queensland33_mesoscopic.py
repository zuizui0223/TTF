#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np


def _cell_key(payload: dict) -> tuple[int, float, float]:
    cell = payload["cell"]
    return (
        int(payload["layout_index"]),
        float(cell["shared_fraction"]),
        float(cell["amplitude"]),
    )


def main() -> int:
    parser = argparse.ArgumentParser(description="Aggregate frozen-layout Queensland mesoscopic calibration.")
    parser.add_argument("--input-dir", type=Path, required=True)
    parser.add_argument("--previous", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--type1-ceiling", type=float, default=0.10)
    parser.add_argument("--power-floor", type=float, default=0.80)
    args = parser.parse_args()

    payloads = []
    for path in sorted(args.input_dir.glob("*.json")):
        payload = json.loads(path.read_text())
        if payload.get("schema") != "ttf_queensland33_mesoscopic_cell_v0.3":
            continue
        payloads.append(payload)
    if not payloads:
        raise RuntimeError("no mesoscopic cell payloads found")

    cells = {_cell_key(payload): payload for payload in payloads}
    expected_layouts = set(range(12))
    layouts = {key[0] for key in cells}
    if layouts != expected_layouts:
        raise RuntimeError(f"expected layouts 0..11, got {sorted(layouts)}")

    zero_amplitudes = (0.5, 1.0, 2.0, 3.0, 4.0)
    partial_shared = (0.25, 0.5, 1.0)
    for layout in sorted(layouts):
        for amplitude in zero_amplitudes:
            if (layout, 0.0, amplitude) not in cells:
                raise RuntimeError(f"missing private cell layout={layout} amplitude={amplitude}")
        for shared in partial_shared:
            if (layout, shared, 2.0) not in cells:
                raise RuntimeError(f"missing shared cell layout={layout} shared={shared}")

    per_layout: dict[str, dict] = {}
    for layout in sorted(layouts):
        zero = [cells[(layout, 0.0, a)]["cell"] for a in zero_amplitudes]
        full = cells[(layout, 1.0, 2.0)]["cell"]
        partial = {
            str(shared): {
                "rejection_rate": float(cells[(layout, shared, 2.0)]["cell"]["rejection_rate"]),
                "mean_statistic": float(cells[(layout, shared, 2.0)]["cell"]["mean_statistic"]),
            }
            for shared in partial_shared
        }
        max_type1 = max(float(cell["rejection_rate"]) for cell in zero)
        power = float(full["rejection_rate"])
        per_layout[str(layout)] = {
            "max_zero_shared_rejection": max_type1,
            "full_shared_A2_power": power,
            "passed": bool(max_type1 <= float(args.type1_ceiling) and power >= float(args.power_floor)),
            "partial_shared_A2": partial,
        }

    worst_type1 = max(v["max_zero_shared_rejection"] for v in per_layout.values())
    worst_power = min(v["full_shared_A2_power"] for v in per_layout.values())

    fraction_summary = {}
    for shared in (0.0, 0.25, 0.5, 1.0):
        relevant = [
            payload["cell"]
            for (layout, sf, amplitude), payload in cells.items()
            if np.isclose(sf, shared) and np.isclose(amplitude, 2.0)
        ]
        fraction_summary[str(shared)] = {
            "mean_rejection_rate_across_layouts": float(np.mean([c["rejection_rate"] for c in relevant])),
            "min_rejection_rate_across_layouts": float(np.min([c["rejection_rate"] for c in relevant])),
            "max_rejection_rate_across_layouts": float(np.max([c["rejection_rate"] for c in relevant])),
            "mean_statistic_across_layouts": float(np.mean([c["mean_statistic"] for c in relevant])),
        }

    previous = json.loads(args.previous.read_text())
    old = previous["qualification"]
    passed = all(v["passed"] for v in per_layout.values())
    status = (
        "PASS_PILOT"
        if passed
        else (
            "FAIL_TYPE1_AND_POWER"
            if worst_type1 > float(args.type1_ceiling) and worst_power < float(args.power_floor)
            else "FAIL_TYPE1"
            if worst_type1 > float(args.type1_ceiling)
            else "FAIL_POWER"
        )
    )

    output = {
        "schema": "ttf_queensland33_mesoscopic_qualification_v0.3",
        "design": {
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
        "response_to_shared_fraction_A2": fraction_summary,
        "per_layout": per_layout,
        "comparison_to_v0.2": {
            "v0.2_worst_A2_power": float(old["worst_layout_full_shared_A2_power"]),
            "v0.3_worst_A2_power": float(worst_power),
            "delta_worst_A2_power": float(worst_power - float(old["worst_layout_full_shared_A2_power"])),
            "v0.2_worst_type1": float(old["worst_layout_max_zero_shared_rejection"]),
            "v0.3_worst_type1": float(worst_type1),
            "delta_worst_type1": float(worst_type1 - float(old["worst_layout_max_zero_shared_rejection"])),
        },
        "interpretation": {
            "estimator": "Mesoscopic boundary-location evidence is learned on the seven shared inter-basin cuts and evaluated by held-out log-score gain over a geometry-conditioned uniform private-cut baseline.",
            "selection": "Primary prior_strength=1.0 and the existing 0.10/0.80 Gate were fixed before this calibration; no empirical genetic outcome or named Queensland barrier was used.",
            "if_pass": "Repeat on an independent world seed and increase world count before any empirical genetic outcome is opened.",
            "if_fail": "Do not tune on empirical outcomes. Diagnose whether failure is chance-commonness calibration or training-field concentration, then qualify any revision on a fresh synthetic seed.",
        },
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(output, indent=2, sort_keys=True) + "\n")
    print(json.dumps(output["qualification"] | output["comparison_to_v0.2"], sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
