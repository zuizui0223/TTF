from __future__ import annotations

import argparse
import json
from pathlib import Path

from ttf.directionality import directional_coupling_test
from ttf.directionality_simulate import simulate_directional_world


def _world_seed(master_seed: int, cell_index: int, replicate: int) -> int:
    return int(master_seed + 10_000_000 * cell_index + replicate)


def run(rule: dict, *, cell_id: str, start: int, count: int) -> dict:
    if rule["status"] != "frozen_before_q4_1_outcomes":
        raise ValueError("Q4.1 rule is not frozen")
    cells = list(rule["mandatory_cells"])
    matches = [(i, cell) for i, cell in enumerate(cells) if cell["id"] == cell_id]
    if len(matches) != 1:
        raise ValueError(f"unknown or duplicate Q4.1 cell: {cell_id}")
    cell_index, cell = matches[0]
    total = int(rule["inference"]["replicates_per_cell"])
    if start < 0 or count < 1 or start + count > total:
        raise ValueError("invalid replicate slice")

    geo = rule["geometry"]
    direction = rule["direction_rule"]
    strict = rule["strict_private_geometry"]
    world_cfg = rule["world_parameters"]
    inf = rule["inference"]

    rows = []
    for replicate in range(start, start + count):
        seed = _world_seed(int(world_cfg["master_seed"]), cell_index, replicate)
        world = simulate_directional_world(
            mode=cell["mode"],
            n_species=int(geo["n_systems"]),
            records_per_species=int(geo["records_per_system"]),
            low_mismatch=float(world_cfg["low_mismatch"]),
            high_mismatch=float(world_cfg["high_mismatch"]),
            noise_sd=float(cell["noise_sd"]),
            transition_width=float(world_cfg["transition_width"]),
            strict_private_arc_halfwidth=float(strict["arc_halfwidth_radians"]),
            seed=seed,
        )
        labels = sorted(x.species for x in world.samples)
        n_train = int(geo["train_systems"])
        n_eval = int(geo["eval_systems"])
        if len(labels) != n_train + n_eval:
            raise RuntimeError("Q4.1 split size does not match frozen system count")
        train = labels[:n_train]
        evaluation = labels[n_train:]

        result = directional_coupling_test(
            world.samples,
            train_species=train,
            eval_species=evaluation,
            k=int(geo["graph_k"]),
            bandwidth=float(geo["bandwidth"]),
            front_edge_fraction=float(direction["front_edge_fraction"]),
            front_window=float(direction["front_window"]),
            sign_consistency_fraction=float(direction["sign_consistency_fraction"]),
            alpha=float(inf["alpha"]),
            n_bootstrap=int(inf["bootstrap_resamples"]),
            seed=seed + 2_000_000,
        )
        delta_values = list(result.species_delta.values())
        rows.append(
            {
                "replicate": int(replicate),
                "world_seed": int(seed),
                "ttf_c_statistic": float(result.coupling.coupling_statistic),
                "ttf_c_p_value": float(result.coupling.coupling_bootstrap.p_value),
                "relation_detected": bool(result.relation_detected),
                "front_center": float(result.front_center),
                "mean_oriented_delta": float(sum(delta_values) / len(delta_values)),
                "breakdown_p_value": float(result.breakdown_bootstrap.p_value),
                "recoupling_p_value": float(result.recoupling_bootstrap.p_value),
                "positive_fraction": float(result.positive_fraction),
                "negative_fraction": float(result.negative_fraction),
                "breakdown_label": bool(result.breakdown_label),
                "recoupling_label": bool(result.recoupling_label),
            }
        )

    return {
        "schema": "ttfm_q4_1_strict_private_batch_v0.1",
        "rule_schema": rule["schema"],
        "cell_id": cell_id,
        "cell_mode": cell["mode"],
        "cell_noise_sd": float(cell["noise_sd"]),
        "cell_roles": list(cell["roles"]),
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
    result = run(rule, cell_id=args.cell_id, start=args.start, count=args.count)
    out = Path(args.output)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps({"cell_id": args.cell_id, "start": args.start, "count": args.count}, sort_keys=True))


if __name__ == "__main__":
    main()
