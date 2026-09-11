from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

from ttf.heterogeneous_inference import heterogeneous_paired_heldout_species_bootstrap_test
from ttf.q5_factor_diagnostic import simulate_factor_world


def _seed(master: int, profile_index: int, response_index: int, replicate: int) -> int:
    return int(master + 10_000_000 * profile_index + 1_000_000 * response_index + replicate)


def run(rule: dict, *, profile: str, response_id: str, start: int, count: int) -> dict:
    if rule["status"] != "frozen_before_development_diagnostic_outcomes":
        raise ValueError("diagnostic rule is not frozen")
    profiles = list(rule["geometry_profiles"])
    if profile not in profiles:
        raise ValueError(f"unknown profile: {profile}")
    responses = list(rule["response_families"])
    matches = [(i, x) for i, x in enumerate(responses) if x["id"] == response_id]
    if len(matches) != 1:
        raise ValueError(f"unknown response id: {response_id}")
    response_index, response = matches[0]
    profile_index = profiles.index(profile)
    diag = rule["diagnostic"]
    total = int(diag["replicates_per_cell"])
    if start < 0 or count < 1 or start + count > total:
        raise ValueError("invalid replicate slice")
    const = rule["geometry_constants"]
    train = [f"sp_{i:03d}" for i in range(20)]
    evaluation = [f"sp_{i:03d}" for i in range(20, 40)]
    rows = []
    for replicate in range(start, start + count):
        seed = _seed(int(diag["master_seed"]), profile_index, response_index, replicate)
        world = simulate_factor_world(
            profile=profile,
            response_mode=response["response_mode"],
            amplitude=float(response["amplitude"]),
            noise_sd=float(response["noise_sd"]),
            transition_width=float(const["transition_width"]),
            seed=seed,
        )
        result = heterogeneous_paired_heldout_species_bootstrap_test(
            world.samples,
            train_species=train,
            eval_species=evaluation,
            graph_fraction=float(const["graph_fraction"]),
            bandwidth=float(const["bandwidth"]),
            prior_strength=float(const["prior_strength"]),
            prior_mean=float(const["prior_mean"]),
            segment_points=int(const["segment_points"]),
            n_bootstrap=int(diag["bootstrap_resamples"]),
            seed=seed + 500_000,
            edge_chunk_size=32,
            train_chunk_size=2048,
        )
        train_n = np.asarray([result.effective_n[s] for s in train], dtype=float)
        eval_n = np.asarray([result.effective_n[s] for s in evaluation], dtype=float)
        train_k = np.asarray([result.graph_k[s] for s in train], dtype=float)
        eval_k = np.asarray([result.graph_k[s] for s in evaluation], dtype=float)
        rows.append({
            "replicate": replicate,
            "seed": seed,
            "ttf_m_statistic": float(result.mismatch_statistic),
            "ttf_m_p_value": float(result.mismatch_bootstrap.p_value),
            "ttf_c_statistic": float(result.coupling_statistic),
            "ttf_c_p_value": float(result.coupling_bootstrap.p_value),
            "train_effective_n_mean": float(train_n.mean()),
            "eval_effective_n_mean": float(eval_n.mean()),
            "train_k_mean": float(train_k.mean()),
            "eval_k_mean": float(eval_k.mean()),
        })
    return {
        "schema": "ttfm_q5_factor_diagnostic_batch_v0.1",
        "rule_schema": rule["schema"],
        "profile": profile,
        "response_id": response_id,
        "response_mode": response["response_mode"],
        "start": start,
        "count": count,
        "rows": rows,
        "claim_firewall": rule["claim_firewall"],
    }


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--rule", required=True)
    p.add_argument("--profile", required=True)
    p.add_argument("--response-id", required=True)
    p.add_argument("--start", type=int, required=True)
    p.add_argument("--count", type=int, required=True)
    p.add_argument("--output", required=True)
    args = p.parse_args()
    rule = json.loads(Path(args.rule).read_text(encoding="utf-8"))
    payload = run(rule, profile=args.profile, response_id=args.response_id, start=args.start, count=args.count)
    out = Path(args.output)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps({"profile": args.profile, "response_id": args.response_id, "start": args.start}))


if __name__ == "__main__":
    main()
