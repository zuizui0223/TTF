from __future__ import annotations

import argparse
import json
from pathlib import Path

from ttf.heterogeneous_simulate import simulate_q5_world
from ttf.q52_transport import q52_geometry_transport_development_test


def _world_seed(master_seed: int, cell_index: int, replicate: int) -> int:
    return int(master_seed + 10_000_000 * cell_index + replicate)


def run(rule: dict, *, cell_id: str, start: int, count: int) -> dict:
    if rule["status"] != "frozen_before_q5_2_transport_development_outcomes":
        raise ValueError("Q5.2 transport development rule is not frozen")
    cells = list(rule["development_cells"])
    matches = [(i, cell) for i, cell in enumerate(cells) if cell["id"] == cell_id]
    if len(matches) != 1:
        raise ValueError(f"unknown or duplicate cell id: {cell_id}")
    cell_index, cell = matches[0]
    inference = rule["development_inference"]
    total = int(inference["worlds_per_cell"])
    if start < 0 or count < 1 or start + count > total:
        raise ValueError("invalid replicate slice")

    fixed = rule["fixed_q5_1_components"]
    world_cfg = rule["world_parameters"]
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
            transition_width=float(world_cfg["transition_width"]),
            seed=seed,
        )
        result = q52_geometry_transport_development_test(
            world.samples,
            train_species=train,
            eval_species=evaluation,
            graph_fraction=float(fixed["graph_fraction"]),
            bandwidth=float(fixed["bandwidth"]),
            prior_strength=float(fixed["prior_strength"]),
            segment_points=int(fixed["segment_points"]),
            n_bootstrap=int(inference["bootstrap_resamples"]),
            seed=seed + 2_000_000,
            edge_chunk_size=32,
            train_chunk_size=2048,
            density_chunk_size=256,
        )
        rows.append(
            {
                "replicate": int(replicate),
                "world_seed": int(seed),
                "modes": {
                    mode: {
                        "statistic": float(mode_result.statistic),
                        "p_value": float(mode_result.bootstrap.p_value),
                        "mean_train_effective_edge_count": float(
                            mode_result.mean_train_effective_edge_count
                        ),
                        "min_train_effective_edge_count": float(
                            mode_result.min_train_effective_edge_count
                        ),
                    }
                    for mode, mode_result in result.modes.items()
                },
            }
        )

    return {
        "schema": "ttfm_q5_2_geometry_transport_development_batch_v0.1",
        "rule_schema": rule["schema"],
        "cell_id": cell_id,
        "response_mode": cell["response_mode"],
        "geometry_profile": cell["geometry_profile"],
        "role": cell["role"],
        "start": int(start),
        "count": int(count),
        "rows": rows,
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
