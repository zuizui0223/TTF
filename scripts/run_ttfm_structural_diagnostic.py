from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

from ttf.core import split_species
from ttf.mismatch import build_coupling_edges, build_mismatch_edges
from ttf.mismatch_simulate import simulate_paired_transition_world
from ttf.transfer import transfer_statistic


def _summary(values: list[float]) -> dict[str, float]:
    x = np.asarray(values, dtype=float)
    return {
        "mean": float(np.mean(x)),
        "median": float(np.median(x)),
        "q05": float(np.quantile(x, 0.05)),
        "q95": float(np.quantile(x, 0.95)),
        "fraction_positive": float(np.mean(x > 0.0)),
    }


def _seed(master: int, mode_index: int, replicate: int) -> int:
    # Fixed arithmetic namespace; no Python hash dependence.
    return int(master + 100_000 * mode_index + replicate)


def run(rule: dict) -> dict:
    modes = list(rule["worlds"].keys())
    out: dict[str, object] = {
        "schema": "ttfm_structural_diagnostic_result_v0.1",
        "status": "development_only_no_formal_pass_fail",
        "rule_status": rule["status"],
        "settings": {
            key: rule[key]
            for key in [
                "n_species",
                "records_per_species",
                "amplitude",
                "noise_sd",
                "transition_width",
                "graph_k",
                "eval_fraction",
                "bandwidth",
                "prior_strength",
                "segment_points",
                "master_seed",
            ]
        },
        "claim_firewall": rule["claim_firewall"],
        "formal_pass_fail": None,
        "modes": {},
    }

    for mode_index, mode in enumerate(modes):
        n_worlds = int(rule["worlds"][mode])
        mismatch_stats: list[float] = []
        coupling_stats: list[float] = []

        for replicate in range(n_worlds):
            world_seed = _seed(int(rule["master_seed"]), mode_index, replicate)
            world = simulate_paired_transition_world(
                mode=mode,
                n_species=int(rule["n_species"]),
                records_per_species=int(rule["records_per_species"]),
                amplitude=float(rule["amplitude"]),
                noise_sd=float(rule["noise_sd"]),
                transition_width=float(rule["transition_width"]),
                seed=world_seed,
            )
            labels = [sample.species for sample in world.samples]
            train_ids, eval_ids = split_species(
                labels,
                eval_fraction=float(rule["eval_fraction"]),
                seed=world_seed + 50_000,
            )
            by_name = {sample.species: sample for sample in world.samples}

            train_m = [
                build_mismatch_edges(by_name[name], k=int(rule["graph_k"]))
                for name in train_ids
            ]
            eval_m = [
                build_mismatch_edges(by_name[name], k=int(rule["graph_k"]))
                for name in eval_ids
            ]
            train_c = [
                build_coupling_edges(by_name[name], k=int(rule["graph_k"]))
                for name in train_ids
            ]
            eval_c = [
                build_coupling_edges(by_name[name], k=int(rule["graph_k"]))
                for name in eval_ids
            ]

            mismatch_stats.append(
                transfer_statistic(
                    train_m,
                    eval_m,
                    bandwidth=float(rule["bandwidth"]),
                    prior_strength=float(rule["prior_strength"]),
                    segment_points=int(rule["segment_points"]),
                ).statistic
            )
            coupling_stats.append(
                transfer_statistic(
                    train_c,
                    eval_c,
                    bandwidth=float(rule["bandwidth"]),
                    prior_strength=float(rule["prior_strength"]),
                    segment_points=int(rule["segment_points"]),
                ).statistic
            )

        out["modes"][mode] = {
            "n_worlds": n_worlds,
            "ttf_m_mismatch_transfer": _summary(mismatch_stats),
            "ttf_c_coupling_transfer": _summary(coupling_stats),
        }

    return out


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--rule", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    rule = json.loads(Path(args.rule).read_text(encoding="utf-8"))
    assert rule["status"] == "frozen_before_structural_diagnostic_outcomes"
    assert rule["formal_pass_fail_forbidden"] is True
    assert rule["candidate_selection_forbidden"] is True
    assert rule["empirical_data_forbidden"] is True

    result = run(rule)
    path = Path(args.output)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
