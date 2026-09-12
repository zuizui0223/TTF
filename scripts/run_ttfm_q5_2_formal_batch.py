from __future__ import annotations

import argparse
import json
from pathlib import Path

from ttf.heterogeneous_simulate import simulate_q5_world
from ttf.q52_formal import q52_selected_transport_test


def _world_seed(master_seed: int, cell_index: int, replicate: int) -> int:
    return int(master_seed + 10_000_000 * cell_index + replicate)


def selected_candidate(protocol: dict, development: dict) -> str:
    if protocol.get("status") != "frozen_before_q5_2_development_aggregate_and_formal_outcomes":
        raise ValueError("Q5.2 conditional formal protocol is not frozen")
    if development.get("schema") != "ttfm_q5_2_geometry_transport_development_result_v0.1":
        raise ValueError("unexpected Q5.2 development result schema")
    candidate = development.get("candidate_selected")
    if candidate is None:
        raise ValueError("Q5.2 development selected no candidate; formal qualification is forbidden")
    allowed = set(protocol["candidate_handoff"]["allowed_values"])
    if candidate not in allowed:
        raise ValueError("development selected a candidate outside the frozen formal protocol")
    checks = development.get("candidate_checks", {})
    if candidate not in checks or checks[candidate].get("development_pass") is not True:
        raise ValueError("selected Q5.2 candidate did not pass the frozen development gate")
    return str(candidate)


def run(
    protocol: dict,
    development: dict,
    *,
    cell_id: str,
    start: int,
    count: int,
) -> dict:
    candidate = selected_candidate(protocol, development)
    cells = list(protocol["formal_cells"])
    matches = [(i, cell) for i, cell in enumerate(cells) if cell["id"] == cell_id]
    if len(matches) != 1:
        raise ValueError(f"unknown or duplicate formal cell id: {cell_id}")
    cell_index, cell = matches[0]
    formal = protocol["formal_worlds"]
    total = int(formal["worlds_per_cell"])
    if start < 0 or count < 1 or start + count > total:
        raise ValueError("invalid formal replicate slice")

    fixed = protocol["fixed_method"]
    train = [f"sp_{i:03d}" for i in range(20)]
    evaluation = [f"sp_{i:03d}" for i in range(20, 40)]
    rows = []
    for replicate in range(start, start + count):
        seed = _world_seed(int(formal["master_seed"]), cell_index, replicate)
        world = simulate_q5_world(
            response_mode=cell["response_mode"],
            geometry_profile=cell["geometry_profile"],
            amplitude=float(cell["amplitude"]),
            noise_sd=float(cell["noise_sd"]),
            transition_width=float(fixed["transition_width"]),
            seed=seed,
        )
        result = q52_selected_transport_test(
            world.samples,
            train_species=train,
            eval_species=evaluation,
            candidate=candidate,
            graph_fraction=float(fixed["graph_fraction"]),
            bandwidth=float(fixed["field_bandwidth"]),
            prior_strength=float(fixed["prior_strength"]),
            segment_points=int(fixed["segment_points"]),
            n_bootstrap=int(formal["bootstrap_resamples"]),
            seed=seed + 2_000_000,
            edge_chunk_size=32,
            train_chunk_size=2048,
            density_chunk_size=256,
        )
        rows.append(
            {
                "replicate": int(replicate),
                "world_seed": int(seed),
                "candidate": candidate,
                "statistic": float(result.statistic),
                "p_value": float(result.bootstrap.p_value),
                "mean_train_effective_edge_count": float(result.mean_train_effective_edge_count),
                "min_train_effective_edge_count": float(result.min_train_effective_edge_count),
            }
        )

    return {
        "schema": "ttfm_q5_2_formal_batch_v0.1",
        "protocol_schema": protocol["schema"],
        "development_result_schema": development["schema"],
        "candidate": candidate,
        "cell_id": cell_id,
        "response_mode": cell["response_mode"],
        "geometry_profile": cell["geometry_profile"],
        "role": cell["role"],
        "start": int(start),
        "count": int(count),
        "rows": rows,
        "claim_firewall": protocol["claim_firewall"],
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--protocol", required=True)
    parser.add_argument("--development-result", required=True)
    parser.add_argument("--cell-id", required=True)
    parser.add_argument("--start", type=int, required=True)
    parser.add_argument("--count", type=int, required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    protocol = json.loads(Path(args.protocol).read_text(encoding="utf-8"))
    development = json.loads(Path(args.development_result).read_text(encoding="utf-8"))
    payload = run(
        protocol,
        development,
        cell_id=args.cell_id,
        start=args.start,
        count=args.count,
    )
    out = Path(args.output)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps({"candidate": payload["candidate"], "cell_id": args.cell_id, "start": args.start, "count": args.count}, sort_keys=True))


if __name__ == "__main__":
    main()
