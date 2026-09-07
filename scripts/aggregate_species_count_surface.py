#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

from ttf.precision import wilson_interval

EXPECTED_TOTAL_SPECIES = (12, 16, 20, 30, 40, 60)
EXPECTED_ARMS = ((0.0, 3.0), (1.0, 2.0))


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Reduce the v0.3 held-out-species count development surface."
    )
    parser.add_argument("--input-dir", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--type1-upper-ceiling", type=float, default=0.10)
    parser.add_argument("--power-lower-floor", type=float, default=0.80)
    args = parser.parse_args()

    paths = sorted(args.input_dir.glob("cell-*.json"))
    expected_count = len(EXPECTED_TOTAL_SPECIES) * len(EXPECTED_ARMS)
    if len(paths) != expected_count:
        raise RuntimeError(f"expected {expected_count} cells, got {len(paths)}")

    rows: dict[int, dict[str, object]] = {
        n: {
            "total_species": n,
            "heldout_species": int(round(n * 0.5)),
            "training_species": n - int(round(n * 0.5)),
        }
        for n in EXPECTED_TOTAL_SPECIES
    }
    seen: set[tuple[int, float, float]] = set()
    frozen_config: dict[str, object] | None = None

    for path in paths:
        payload = json.loads(path.read_text())
        if payload.get("schema") != "ttf_calibration_cell_v0.1":
            raise RuntimeError(f"unexpected cell schema in {path}")
        config = payload["config"]
        cell = payload["cell"]
        total = int(config["species"])
        shared = float(cell["shared_fraction"])
        amplitude = float(cell["amplitude"])
        key = (total, shared, amplitude)
        if key in seen:
            raise RuntimeError(f"duplicate surface cell {key}")
        seen.add(key)
        if total not in rows or (shared, amplitude) not in EXPECTED_ARMS:
            raise RuntimeError(f"unexpected surface cell {key}")

        invariant = {
            "replicates": int(config["replicates"]),
            "resamples": int(config["permutations"]),
            "records_per_species": int(config["records_per_species"]),
            "bandwidth": float(config["bandwidth"]),
            "alpha": float(config["alpha"]),
            "seed": int(config["seed"]),
            "inference": str(config["inference"]),
        }
        if frozen_config is None:
            frozen_config = invariant
        elif invariant != frozen_config:
            raise RuntimeError("species-count surface config drift")

        trials = int(cell["n_replicates"])
        successes_raw = float(cell["rejection_rate"]) * trials
        successes = int(round(successes_raw))
        if not np.isclose(successes_raw, successes, atol=1e-8, rtol=0.0):
            raise RuntimeError("rejection rate does not map to integer successes")
        interval = wilson_interval(successes, trials).to_dict()
        entry = {
            "shared_fraction": shared,
            "amplitude": amplitude,
            "mean_statistic": float(cell["mean_statistic"]),
            "median_p_value": float(cell["median_p_value"]),
            **interval,
        }
        if shared == 0.0:
            rows[total]["zero_shared_amplitude_3"] = entry
        else:
            rows[total]["full_shared_amplitude_2"] = entry

    expected = {
        (n, shared, amplitude)
        for n in EXPECTED_TOTAL_SPECIES
        for shared, amplitude in EXPECTED_ARMS
    }
    if seen != expected:
        raise RuntimeError("species-count surface cell set is incomplete")
    assert frozen_config is not None
    if frozen_config != {
        "replicates": 500,
        "resamples": 1999,
        "records_per_species": 60,
        "bandwidth": 0.2,
        "alpha": 0.05,
        "seed": 20260908,
        "inference": "heldout_species_bootstrap",
    }:
        raise RuntimeError(f"surface design drift: {frozen_config}")

    ordered: list[dict[str, object]] = []
    for total in EXPECTED_TOTAL_SPECIES:
        row = rows[total]
        zero = row["zero_shared_amplitude_3"]
        positive = row["full_shared_amplitude_2"]
        type1_pass = bool(float(zero["high"]) <= args.type1_upper_ceiling)
        power_pass = bool(float(positive["low"]) >= args.power_lower_floor)
        row["type1_pass"] = type1_pass
        row["power_pass"] = power_pass
        row["joint_pass"] = bool(type1_pass and power_pass)
        ordered.append(row)

    # Prospective stability rule: the floor is the smallest tested held-out count
    # whose cell passes both gates AND every larger tested held-out count also passes.
    recommended: int | None = None
    for i, row in enumerate(ordered):
        if all(bool(later["joint_pass"]) for later in ordered[i:]):
            recommended = int(row["heldout_species"])
            break

    payload = {
        "schema": "ttf_v03_species_count_surface_v0.1",
        "status": "development_surface_complete",
        "trigger": (
            "Gate-I-A v0.1 used 8 held-out species and failed both high-precision "
            "type-I and power gates; v0.2 idealized qualification used 20 held-out species."
        ),
        "estimator_changed": False,
        "development_only": True,
        "claim_ready": False,
        "design": {
            **frozen_config,
            "total_species_grid": list(EXPECTED_TOTAL_SPECIES),
            "heldout_species_grid": [int(round(n * 0.5)) for n in EXPECTED_TOTAL_SPECIES],
            "eval_fraction": 0.5,
            "zero_shared_stress_amplitude": 3.0,
            "positive_control_amplitude": 2.0,
            "type1_upper_ceiling": args.type1_upper_ceiling,
            "power_lower_floor": args.power_lower_floor,
            "interval": "two-sided Wilson score 95%",
            "floor_selection_rule": (
                "smallest tested held-out species count that passes both gates and for "
                "which every larger tested held-out count also passes both gates"
            ),
        },
        "rows": ordered,
        "recommended_min_claim_heldout_species": recommended,
        "recommended_min_claim_training_species": recommended,
        "successor_requirement": (
            "Any non-null recommended floor is a development result only. It must be "
            "frozen before acquiring/selecting the fresh external I-A2 confirmation "
            "panel, and I-A v0.1 cannot be reused as confirmatory evidence."
        ),
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    print(json.dumps({
        "recommended_min_claim_heldout_species": recommended,
        "rows": [
            {
                "heldout": row["heldout_species"],
                "type1_pass": row["type1_pass"],
                "power_pass": row["power_pass"],
            }
            for row in ordered
        ],
    }, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
