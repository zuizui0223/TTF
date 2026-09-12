from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

from ttf.heterogeneous_simulate import simulate_q5_world
from ttf.support_overlap_diagnostic import paired_support_overlap_diagnostic


def _world_seed(master_seed: int, cell_index: int, replicate: int) -> int:
    return int(master_seed + 10_000_000 * cell_index + replicate)


def _scalar(array: np.ndarray) -> float:
    values = np.asarray(array, dtype=float)
    if values.shape != (1,) or not np.isfinite(values[0]):
        raise ValueError("diagnostic batch expects one finite world per scorer call")
    return float(values[0])


def run(rule: dict, *, cell_id: str, start: int, count: int) -> dict:
    if rule.get("status") != "frozen_before_support_overlap_diagnostic_outcomes":
        raise ValueError("support-overlap diagnostic rule is not frozen")
    if rule.get("not_a_candidate_family") is not True:
        raise ValueError("support-overlap diagnostic must remain non-candidate")

    cells = list(rule["cells"])
    matches = [(i, cell) for i, cell in enumerate(cells) if cell["id"] == cell_id]
    if len(matches) != 1:
        raise ValueError(f"unknown or duplicate cell id: {cell_id}")
    cell_index, cell = matches[0]

    world_cfg = rule["diagnostic_worlds"]
    total = int(world_cfg["worlds_per_cell"])
    if start < 0 or count < 1 or start + count > total:
        raise ValueError("invalid replicate slice")

    fixed = rule["fixed_estimators"]
    train = [f"sp_{i:03d}" for i in range(20)]
    evaluation = [f"sp_{i:03d}" for i in range(20, 40)]
    rows = []

    for replicate in range(start, start + count):
        seed = _world_seed(int(world_cfg["master_seed"]), cell_index, replicate)
        world = simulate_q5_world(
            response_mode=cell["response_mode"],
            geometry_profile=cell["geometry_profile"],
            amplitude=float(cell["amplitude"]),
            noise_sd=float(cell["noise_sd"]),
            transition_width=0.2,
            seed=seed,
        )
        result = paired_support_overlap_diagnostic(
            world.samples,
            train_species=train,
            eval_species=evaluation,
            graph_fraction=float(fixed["graph_fraction"]),
            bandwidth=float(fixed["bandwidth"]),
            prior_strength=float(fixed["prior_strength"]),
            segment_points=int(fixed["segment_points"]),
            edge_chunk_size=32,
            train_chunk_size=2048,
        )

        # Geometry-only support summaries are shared exactly between M and C.
        support = result.coupling
        eval_systems = {}
        for name in evaluation:
            eval_systems[name] = {
                "mismatch": {
                    "raw_score": _scalar(result.mismatch.species_raw_scores[name]),
                    "opportunity_partial_score": _scalar(
                        result.mismatch.species_partial_scores[name]
                    ),
                    "prediction_opportunity_rho": _scalar(
                        result.mismatch.species_prediction_opportunity_rho[name]
                    ),
                    "target_opportunity_rho": _scalar(
                        result.mismatch.species_target_opportunity_rho[name]
                    ),
                    "alignment_product": _scalar(
                        result.mismatch.species_alignment_product[name]
                    ),
                },
                "coupling": {
                    "raw_score": _scalar(result.coupling.species_raw_scores[name]),
                    "opportunity_partial_score": _scalar(
                        result.coupling.species_partial_scores[name]
                    ),
                    "prediction_opportunity_rho": _scalar(
                        result.coupling.species_prediction_opportunity_rho[name]
                    ),
                    "target_opportunity_rho": _scalar(
                        result.coupling.species_target_opportunity_rho[name]
                    ),
                    "alignment_product": _scalar(
                        result.coupling.species_alignment_product[name]
                    ),
                },
                "support": {
                    "mean_opportunity": float(support.species_mean_opportunity[name]),
                    "mean_prior_fraction": float(
                        support.species_mean_prior_fraction[name]
                    ),
                    "low_support_fraction": float(
                        support.species_low_support_fraction[name]
                    ),
                },
            }

        rows.append(
            {
                "replicate": int(replicate),
                "world_seed": int(seed),
                "world": {
                    "raw_M": _scalar(result.mismatch.raw_statistics),
                    "partial_M": _scalar(
                        result.mismatch.opportunity_partial_statistics
                    ),
                    "alignment_M": _scalar(result.mismatch.alignment_statistics),
                    "raw_C": _scalar(result.coupling.raw_statistics),
                    "partial_C": _scalar(
                        result.coupling.opportunity_partial_statistics
                    ),
                    "alignment_C": _scalar(result.coupling.alignment_statistics),
                    "mean_opportunity": float(support.mean_opportunity),
                    "mean_prior_fraction": float(support.mean_prior_fraction),
                    "low_support_fraction": float(support.low_support_fraction),
                },
                "eval_systems": eval_systems,
                "graph_k": {str(k): int(v) for k, v in result.graph_k.items()},
                "effective_n": {
                    str(k): int(v) for k, v in result.effective_n.items()
                },
            }
        )

    return {
        "schema": "ttfc_support_overlap_diagnostic_batch_v0.1",
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
    print(
        json.dumps(
            {"cell_id": args.cell_id, "start": args.start, "count": args.count},
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()
