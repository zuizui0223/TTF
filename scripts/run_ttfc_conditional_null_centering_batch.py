from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

from ttf.conditional_null_centering import (
    identical_uniform_circle_samples,
    q51_private_relation_phase_redraw_scores,
)
from ttf.heterogeneous_simulate import simulate_q5_world


def _geometry_seed(master: int, cell_index: int, replicate: int) -> int:
    return int(master + 10_000_000 * cell_index + replicate)


def _phase_seed(master: int, cell_index: int, replicate: int) -> int:
    return int(master + 10_000_000 * cell_index + replicate)


def _samples_for_cell(cell: dict, seed: int):
    mode = str(cell["geometry_mode"])
    if mode == "identical_uniform60":
        return identical_uniform_circle_samples(seed=int(seed), n_systems=40, n_points=60)
    if mode in {"q5_matched", "q5_shifted"}:
        profile = "matched" if mode == "q5_matched" else "shifted"
        # Q5 realizes all geometry before phases/responses. Those generated
        # responses are discarded; only coordinates are passed downstream.
        world = simulate_q5_world(
            response_mode="null_component",
            geometry_profile=profile,
            amplitude=2.0,
            noise_sd=0.8,
            transition_width=0.2,
            seed=int(seed),
        )
        return world.samples
    raise ValueError(f"unknown geometry_mode: {mode}")


def run(rule: dict, *, cell_id: str, start: int, count: int) -> dict:
    if rule.get("status") != "frozen_before_conditional_null_centering_outcomes":
        raise ValueError("conditional-null centering rule is not frozen")
    if rule.get("not_a_candidate_family") is not True:
        raise ValueError("conditional-null centering must remain non-candidate")

    cells = list(rule["cells"])
    matches = [(i, cell) for i, cell in enumerate(cells) if cell["id"] == cell_id]
    if len(matches) != 1:
        raise ValueError(f"unknown or duplicate cell id: {cell_id}")
    cell_index, cell = matches[0]

    cfg = rule["diagnostic_worlds"]
    total = int(cfg["geometry_worlds_per_cell"])
    if start < 0 or count < 1 or start + count > total:
        raise ValueError("invalid replicate slice")

    fixed = rule["fixed_estimator"]
    train = [f"sp_{i:03d}" for i in range(20)]
    evaluation = [f"sp_{i:03d}" for i in range(20, 40)]
    n_phase = int(cfg["private_phase_draws_per_geometry"])
    rows = []

    for replicate in range(start, start + count):
        geometry_seed = _geometry_seed(
            int(cfg["geometry_master_seed"]), cell_index, replicate
        )
        phase_seed = _phase_seed(
            int(cfg["phase_master_seed"]), cell_index, replicate
        )
        samples = _samples_for_cell(cell, geometry_seed)
        result = q51_private_relation_phase_redraw_scores(
            samples,
            train_species=train,
            eval_species=evaluation,
            n_phase_draws=n_phase,
            phase_seed=phase_seed,
            graph_fraction=float(fixed["graph_fraction"]),
            bandwidth=float(fixed["bandwidth"]),
            prior_strength=float(fixed["prior_strength"]),
            segment_points=int(fixed["segment_points"]),
            edge_chunk_size=32,
            train_chunk_size=2048,
        )
        scores = np.asarray(result.statistics, dtype=float)
        if scores.shape != (n_phase,) or not np.isfinite(scores).all():
            raise RuntimeError("conditional-null score width drift")
        sd = float(np.std(scores, ddof=1))
        rows.append(
            {
                "replicate": int(replicate),
                "geometry_seed": int(geometry_seed),
                "phase_seed": int(phase_seed),
                "conditional_mean_C": float(np.mean(scores)),
                "within_geometry_phase_sd_C": sd,
                "conditional_mean_mc_se_C": float(sd / np.sqrt(len(scores))),
                "phase_fraction_positive_C": float(np.mean(scores > 0.0)),
                "phase_q05_C": float(np.quantile(scores, 0.05)),
                "phase_q50_C": float(np.quantile(scores, 0.50)),
                "phase_q95_C": float(np.quantile(scores, 0.95)),
                "graph_k": {str(k): int(v) for k, v in result.graph_k.items()},
                "effective_n": {
                    str(k): int(v) for k, v in result.effective_n.items()
                },
            }
        )

    return {
        "schema": "ttfc_conditional_null_centering_batch_v0.1",
        "rule_schema": rule["schema"],
        "cell_id": str(cell_id),
        "geometry_mode": str(cell["geometry_mode"]),
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
