from __future__ import annotations

import argparse
import json
from pathlib import Path

from ttf.core import split_species
from ttf.mismatch_inference import paired_heldout_species_bootstrap_test
from ttf.mismatch_simulate import simulate_paired_transition_world


def _world_seed(master_seed: int, cell_index: int, replicate: int) -> int:
    return int(master_seed + 10_000_000 * cell_index + replicate)


def run(rule: dict, *, cell_id: str, start: int, count: int) -> dict:
    if rule["status"] != "frozen_before_formal_qualification_outcomes":
        raise ValueError("formal qualification rule is not frozen")
    cells = list(rule["mandatory_cells"])
    matches = [(i, cell) for i, cell in enumerate(cells) if cell["id"] == cell_id]
    if len(matches) != 1:
        raise ValueError(f"unknown or duplicate cell id: {cell_id}")
    cell_index, cell = matches[0]

    total = int(rule["inference"]["replicates_per_cell"])
    if start < 0 or count < 1 or start + count > total:
        raise ValueError("invalid replicate slice")

    geo = rule["geometry"]
    world_cfg = rule["world_parameters"]
    inf = rule["inference"]
    rows = []
    for replicate in range(start, start + count):
        seed = _world_seed(int(world_cfg["master_seed"]), cell_index, replicate)
        world = simulate_paired_transition_world(
            mode=cell["mode"],
            n_species=int(geo["n_systems"]),
            records_per_species=int(geo["records_per_system"]),
            amplitude=float(cell["amplitude"]),
            noise_sd=float(cell["noise_sd"]),
            transition_width=float(world_cfg["transition_width"]),
            seed=seed,
        )
        labels = [sample.species for sample in world.samples]
        train, evaluation = split_species(
            labels,
            eval_fraction=float(geo["eval_fraction"]),
            seed=seed + 1_000_000,
        )
        result = paired_heldout_species_bootstrap_test(
            world.samples,
            train_species=train,
            eval_species=evaluation,
            k=int(geo["graph_k"]),
            bandwidth=float(geo["bandwidth"]),
            prior_strength=float(geo["prior_strength"]),
            prior_mean=float(geo["prior_mean"]),
            segment_points=int(geo["segment_points"]),
            n_bootstrap=int(inf["bootstrap_resamples"]),
            seed=seed + 2_000_000,
            edge_chunk_size=32,
            train_chunk_size=2048,
        )
        rows.append(
            {
                "replicate": int(replicate),
                "world_seed": int(seed),
                "ttf_m_statistic": float(result.mismatch_statistic),
                "ttf_m_p_value": float(result.mismatch_bootstrap.p_value),
                "ttf_c_statistic": float(result.coupling_statistic),
                "ttf_c_p_value": float(result.coupling_bootstrap.p_value),
            }
        )

    return {
        "schema": "ttfm_idealized_qualification_batch_v0.1",
        "rule_schema": rule["schema"],
        "cell_id": cell_id,
        "cell_mode": cell["mode"],
        "cell_amplitude": float(cell["amplitude"]),
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
