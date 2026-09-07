#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
from pathlib import Path

N_LAYOUTS = 12


def main() -> int:
    parser = argparse.ArgumentParser(description="Aggregate Queensland true-boundary oracle diagnostics.")
    parser.add_argument("--input-dir", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    cells: dict[int, dict[float, dict]] = {i: {} for i in range(N_LAYOUTS)}
    for path in sorted(args.input_dir.rglob("*.json")):
        payload = json.loads(path.read_text())
        if payload.get("schema") != "ttf_queensland33_oracle_diagnostic_cell_v0.1":
            continue
        layout = int(payload["layout_index"])
        sf = float(payload["shared_fraction"])
        if sf in cells[layout]:
            raise RuntimeError(f"duplicate oracle cell layout={layout}, shared={sf}")
        cells[layout][sf] = payload
    for layout in range(N_LAYOUTS):
        if set(cells[layout]) != {0.0, 1.0}:
            raise RuntimeError(f"layout {layout} oracle grid drift")

    per_layout = {}
    worst = {
        "oracle_all_private_rejection": 0.0,
        "oracle_crossing_private_rejection": 0.0,
        "oracle_all_shared_power": 1.0,
        "oracle_crossing_shared_power": 1.0,
    }
    for layout in range(N_LAYOUTS):
        null = cells[layout][0.0]
        shared = cells[layout][1.0]
        all_null = float(null["oracle_all"]["rejection_rate"])
        cross_null = float(null["oracle_crossing"]["rejection_rate"])
        all_power = float(shared["oracle_all"]["rejection_rate"])
        cross_power = float(shared["oracle_crossing"]["rejection_rate"])
        worst["oracle_all_private_rejection"] = max(worst["oracle_all_private_rejection"], all_null)
        worst["oracle_crossing_private_rejection"] = max(worst["oracle_crossing_private_rejection"], cross_null)
        worst["oracle_all_shared_power"] = min(worst["oracle_all_shared_power"], all_power)
        worst["oracle_crossing_shared_power"] = min(worst["oracle_crossing_shared_power"], cross_power)
        per_layout[str(layout)] = {
            "private": {
                "oracle_all_rejection": all_null,
                "oracle_crossing_rejection": cross_null,
            },
            "shared": {
                "oracle_all_power": all_power,
                "oracle_crossing_power": cross_power,
                "oracle_all_mean_score": shared["oracle_all"]["mean_species_score"],
                "oracle_crossing_mean_score": shared["oracle_crossing"]["mean_species_score"],
                "mean_crossing_eval_species": shared["oracle_crossing"]["mean_n_species"],
                "min_crossing_eval_species": shared["oracle_crossing"]["min_n_species"],
            },
        }

    interpretation = {}
    if worst["oracle_all_shared_power"] >= 0.80:
        interpretation["field_learning"] = "The held-out scorer can recover the true transition when handed an oracle boundary; training-field estimation is the main power bottleneck."
    elif worst["oracle_crossing_shared_power"] >= 0.80:
        interpretation["field_learning"] = "Oracle scoring is adequate only after geometry-only observability restriction; both non-observability and training-field estimation matter."
    else:
        interpretation["field_learning"] = "Even an oracle boundary fails the all-layout power gate; sparse held-out edge scoring/inference is itself a structural bottleneck."
    if max(worst["oracle_all_private_rejection"], worst["oracle_crossing_private_rejection"]) > 0.10:
        interpretation["null"] = "A zero-centered species-bootstrap null remains anti-conservative for some coarse discrete layouts even when the oracle cut is independent of private traits; geometry-conditioned chance commonness must be modeled explicitly."
    else:
        interpretation["null"] = "Oracle private-null point rejection remains within the prospective ceiling; excess false positives arise primarily in learned-field transfer."

    payload = {
        "schema": "ttf_queensland33_oracle_diagnostic_v0.1",
        "design": {
            "diagnostic_only": True,
            "n_layouts": N_LAYOUTS,
            "worlds_per_cell_per_layout": 100,
            "bootstrap_resamples": 1999,
            "noise_sd": 0.0,
            "genetic_outcomes_used": False,
            "named_Mary_Brisbane_boundary_used": False,
        },
        "worst_layout_summary": worst,
        "interpretation": interpretation,
        "per_layout": per_layout,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"worst": worst, "interpretation": interpretation}, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
