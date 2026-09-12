from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

from ttf.heterogeneous_simulate import simulate_q5_world
from ttf.training_leverage_diagnostic import paired_training_leverage_diagnostic


def _world_seed(master_seed: int, cell_index: int, replicate: int) -> int:
    return int(master_seed + 10_000_000 * cell_index + replicate)


def _scalar(array: np.ndarray) -> float:
    values = np.asarray(array, dtype=float)
    if values.shape != (1,) or not np.isfinite(values[0]):
        raise ValueError("training-leverage batch expects one finite world per scorer call")
    return float(values[0])


def run(rule: dict, *, cell_id: str, start: int, count: int) -> dict:
    if rule.get("status") != "frozen_before_training_leverage_diagnostic_outcomes":
        raise ValueError("training-leverage diagnostic rule is not frozen")
    if rule.get("not_a_candidate_family") is not True:
        raise ValueError("training-leverage diagnostic must remain non-candidate")

    cells = list(rule["cells"])
    matches = [(i, cell) for i, cell in enumerate(cells) if cell["id"] == cell_id]
    if len(matches) != 1:
        raise ValueError(f"unknown or duplicate cell id: {cell_id}")
    cell_index, cell = matches[0]

    cfg = rule["diagnostic_worlds"]
    total = int(cfg["worlds_per_cell"])
    if start < 0 or count < 1 or start + count > total:
        raise ValueError("invalid replicate slice")

    fixed = rule["fixed_estimators"]
    train = [f"sp_{i:03d}" for i in range(20)]
    evaluation = [f"sp_{i:03d}" for i in range(20, 40)]
    rows = []

    for replicate in range(start, start + count):
        seed = _world_seed(int(cfg["master_seed"]), cell_index, replicate)
        world = simulate_q5_world(
            response_mode=cell["response_mode"],
            geometry_profile=cell["geometry_profile"],
            amplitude=float(cell["amplitude"]),
            noise_sd=float(cell["noise_sd"]),
            transition_width=0.2,
            seed=seed,
        )
        result = paired_training_leverage_diagnostic(
            world.samples,
            train_species=train,
            eval_species=evaluation,
            graph_fraction=float(fixed["graph_fraction"]),
            bandwidth=float(fixed["bandwidth"]),
            prior_strength=float(fixed["prior_strength"]),
            segment_points=int(fixed["segment_points"]),
        )

        loo_by_training_system = {
            name: {
                "M": _scalar(result.mismatch.loo_statistics_by_training_system[name]),
                "C": _scalar(result.coupling.loo_statistics_by_training_system[name]),
            }
            for name in train
        }

        rows.append(
            {
                "replicate": int(replicate),
                "world_seed": int(seed),
                "world": {
                    "raw_M": _scalar(result.mismatch.raw_statistics),
                    "raw_C": _scalar(result.coupling.raw_statistics),
                    "loo_score_sd_M": _scalar(result.mismatch.loo_score_sd),
                    "loo_score_sd_C": _scalar(result.coupling.loo_score_sd),
                    "loo_score_mean_M": _scalar(result.mismatch.loo_score_mean),
                    "loo_score_mean_C": _scalar(result.coupling.loo_score_mean),
                    "loo_score_min_M": _scalar(result.mismatch.loo_score_min),
                    "loo_score_min_C": _scalar(result.coupling.loo_score_min),
                    "mean_abs_loo_score_shift_M": _scalar(
                        result.mismatch.mean_abs_loo_score_shift
                    ),
                    "mean_abs_loo_score_shift_C": _scalar(
                        result.coupling.mean_abs_loo_score_shift
                    ),
                    "mean_prediction_loo_sd_M": _scalar(
                        result.mismatch.mean_prediction_loo_sd
                    ),
                    "mean_prediction_loo_sd_C": _scalar(
                        result.coupling.mean_prediction_loo_sd
                    ),
                    "mean_effective_training_system_count": float(
                        result.mean_effective_training_system_count
                    ),
                    "mean_dominant_training_system_share": float(
                        result.mean_dominant_training_system_share
                    ),
                },
                "loo_by_training_system": loo_by_training_system,
                "graph_k": {str(k): int(v) for k, v in result.graph_k.items()},
                "effective_n": {str(k): int(v) for k, v in result.effective_n.items()},
            }
        )

    return {
        "schema": "ttfc_training_leverage_diagnostic_batch_v0.1",
        "rule_schema": rule["schema"],
        "cell_id": cell_id,
        "response_mode": cell["response_mode"],
        "geometry_profile": cell["geometry_profile"],
        "start": int(start),
        "count": int(count),
        "rows": rows,
        "candidate_selected": None,
        "formal_pass_fail": None,
        "claim_firewall": rule["claim_firewall"],
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--rule", required=True)
    parser.add_argument("--cell-id", required=True)
    parser.add_argument("--start", type=int, required=True)
    parser.add_argument("--count", type=int, required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    rule = json.loads(Path(args.rule).read_text(encoding="utf-8"))
    payload = run(rule, cell_id=args.cell_id, start=args.start, count=args.count)
    out = Path(args.output)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps({"cell_id": args.cell_id, "start": args.start, "count": args.count}, sort_keys=True))


if __name__ == "__main__":
    main()
