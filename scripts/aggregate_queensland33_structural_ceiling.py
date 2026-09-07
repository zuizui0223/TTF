#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
from pathlib import Path

N_LAYOUTS = 12


def main() -> int:
    parser = argparse.ArgumentParser(description="Aggregate the Queensland noise-free structural ceiling.")
    parser.add_argument("--input-dir", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    cells: dict[int, dict[float, dict]] = {i: {} for i in range(N_LAYOUTS)}
    coverage: dict[int, dict | None] = {}
    for path in sorted(args.input_dir.rglob("*.json")):
        payload = json.loads(path.read_text())
        if payload.get("schema") != "ttf_queensland33_structural_ceiling_cell_v0.1":
            continue
        layout = int(payload["layout_index"])
        sf = float(payload["cell"]["shared_fraction"])
        if sf in cells[layout]:
            raise RuntimeError(f"duplicate ceiling cell layout={layout}, shared={sf}")
        cells[layout][sf] = payload["cell"]
        if sf == 1.0:
            coverage[layout] = payload.get("coverage")

    for layout in range(N_LAYOUTS):
        if set(cells[layout]) != {0.0, 1.0}:
            raise RuntimeError(f"layout {layout} ceiling grid drift: {sorted(cells[layout])}")

    per_layout = {}
    worst_private = 0.0
    worst_power = 1.0
    min_eval_crossing_fraction = 1.0
    for layout in range(N_LAYOUTS):
        private = float(cells[layout][0.0]["rejection_rate"])
        power = float(cells[layout][1.0]["rejection_rate"])
        cov = coverage.get(layout)
        crossing_fraction = None if cov is None else float(cov["mean_eval_crossing_fraction"])
        worst_private = max(worst_private, private)
        worst_power = min(worst_power, power)
        if crossing_fraction is not None:
            min_eval_crossing_fraction = min(min_eval_crossing_fraction, crossing_fraction)
        per_layout[str(layout)] = {
            "private_noise_free_rejection": private,
            "full_shared_noise_free_power": power,
            "coverage": cov,
        }

    payload = {
        "schema": "ttf_queensland33_structural_ceiling_v0.1",
        "design": {
            "exact_frozen_layout_envelope": True,
            "n_layouts": N_LAYOUTS,
            "noise_sd": 0.0,
            "amplitude": 1.0,
            "amplitude_interpretation": "scale cancels under within-species ranks at zero noise; this estimates a practical infinite-SNR ceiling",
            "worlds_per_cell_per_layout": 100,
            "bootstrap_resamples": 1999,
            "same_estimator": True,
            "same_10_train_11_eval_split": True,
            "same_k": 2,
            "same_transition_width": 0.2,
            "same_min_shared_crossing_species": 14,
            "genetic_outcomes_used": False,
            "named_Mary_Brisbane_boundary_used": False,
        },
        "summary": {
            "worst_layout_private_noise_free_rejection": worst_private,
            "worst_layout_full_shared_noise_free_power": worst_power,
            "minimum_mean_eval_crossing_fraction_across_layouts": min_eval_crossing_fraction,
            "structurally_qualified_at_0.80": bool(worst_power >= 0.80 and worst_private <= 0.10),
        },
        "interpretation": {
            "if_power_below_0.80": "The frozen sampling topology itself limits held-out detectability even at zero measurement noise; do not rescue the benchmark by increasing trait amplitude.",
            "if_eval_crossing_low": "Evaluate a prospectively defined geometry-only observability rule as a separate estimator version; species that do not span a learned boundary cannot supply evidence about turnover at that boundary.",
        },
        "per_layout": per_layout,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    print(json.dumps(payload["summary"], sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
