from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

from ttf.conditional_null_centering import identical_uniform_circle_samples
from ttf.heterogeneous_simulate import simulate_q5_world
from ttf.phase_propensity_diagnostic import phase_propensity_phase_redraw_scores


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
    if rule.get("status") != "frozen_before_phase_propensity_outcomes":
        raise ValueError("phase-propensity rule is not frozen")
    if rule["terminal_rule"].get("candidate_selected") is not None:
        raise ValueError("phase-propensity diagnostic must remain non-candidate")
    if rule["terminal_rule"].get("formal_pass_fail") is not None:
        raise ValueError("phase-propensity diagnostic must remain non-qualifying")

    cells = list(rule["fresh_worlds"]["cells"])
    matches = [(i, cell) for i, cell in enumerate(cells) if cell["id"] == cell_id]
    if len(matches) != 1:
        raise ValueError(f"unknown or duplicate cell id: {cell_id}")
    cell_index, cell = matches[0]

    cfg = rule["fresh_worlds"]
    total = int(cfg["worlds_per_cell"])
    if start < 0 or count < 1 or start + count > total:
        raise ValueError("invalid replicate slice")

    fixed = rule["estimators"]["baseline_unscaled_q51"]
    train = [f"sp_{i:03d}" for i in range(20)]
    evaluation = [f"sp_{i:03d}" for i in range(20, 40)]
    n_phase = int(cfg["phase_redraws_per_geometry"])
    rows = []

    for replicate in range(start, start + count):
        geometry_seed = _geometry_seed(int(cfg["geometry_master_seed"]), cell_index, replicate)
        phase_seed = _phase_seed(int(cfg["phase_master_seed"]), cell_index, replicate)
        samples = _samples_for_cell(cell, geometry_seed)
        result = phase_propensity_phase_redraw_scores(
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
        baseline = np.asarray(result.baseline_statistics, dtype=float)
        centered = np.asarray(result.propensity_centered_statistics, dtype=float)
        if baseline.shape != (n_phase,) or centered.shape != (n_phase,):
            raise RuntimeError("phase-propensity score width drift")
        if not np.isfinite(baseline).all() or not np.isfinite(centered).all():
            raise RuntimeError("phase-propensity scores must be finite")
        expected_rms = np.asarray(
            [
                float(np.sqrt(np.mean(np.square(result.expected_training_residual[name]))))
                for name in train
            ],
            dtype=float,
        )
        rows.append(
            {
                "replicate": int(replicate),
                "geometry_seed": int(geometry_seed),
                "phase_seed": int(phase_seed),
                "conditional_mean_baseline_C": float(np.mean(baseline)),
                "conditional_mean_propensity_centered_C": float(np.mean(centered)),
                "conditional_mean_baseline_minus_centered_C": float(np.mean(baseline - centered)),
                "within_geometry_phase_sd_baseline_C": float(np.std(baseline, ddof=1)),
                "within_geometry_phase_sd_propensity_centered_C": float(np.std(centered, ddof=1)),
                "phase_fraction_positive_baseline_C": float(np.mean(baseline > 0.0)),
                "phase_fraction_positive_propensity_centered_C": float(np.mean(centered > 0.0)),
                "mean_training_expected_residual_rms": float(np.mean(expected_rms)),
                "graph_k": {str(k): int(v) for k, v in result.graph_k.items()},
                "effective_n": {str(k): int(v) for k, v in result.effective_n.items()},
            }
        )

    return {
        "schema": "ttfc_phase_propensity_diagnostic_batch_v0.1",
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
