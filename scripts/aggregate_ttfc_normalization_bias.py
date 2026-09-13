from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

WORLD_KEYS = (
    "conditional_mean_normalized_C",
    "conditional_mean_unscaled_C",
    "conditional_mean_normalized_minus_unscaled_C",
    "within_geometry_phase_sd_normalized_C",
    "within_geometry_phase_sd_unscaled_C",
    "phase_fraction_positive_normalized_C",
    "phase_fraction_positive_unscaled_C",
)


def _interval(samples: np.ndarray, estimate: float) -> dict:
    values = np.asarray(samples, dtype=float)
    if values.ndim != 1 or len(values) < 1 or not np.isfinite(values).all():
        raise ValueError("bootstrap samples must be finite and non-empty")
    return {
        "estimate": float(estimate),
        "low95": float(np.quantile(values, 0.025)),
        "high95": float(np.quantile(values, 0.975)),
    }


def _bootstrap_mean(values: np.ndarray, *, rng, n_bootstrap: int) -> dict:
    x = np.asarray(values, dtype=float)
    if x.ndim != 1 or len(x) < 2 or not np.isfinite(x).all():
        raise ValueError("mean bootstrap requires finite 1D values")
    draws = np.empty(int(n_bootstrap), dtype=float)
    for i in range(int(n_bootstrap)):
        idx = rng.integers(0, len(x), size=len(x))
        draws[i] = float(np.mean(x[idx]))
    return _interval(draws, float(np.mean(x)))


def _collect(rule: dict, paths: list[Path]) -> dict[str, dict[str, np.ndarray]]:
    cfg = rule["fresh_worlds"]
    expected_n = int(cfg["worlds_per_cell"])
    cell_rule = {cell["id"]: cell for cell in cfg["cells"]}
    rows_by_cell = {cell_id: {} for cell_id in cell_rule}

    for path in paths:
        payload = json.loads(path.read_text(encoding="utf-8"))
        if payload.get("schema") != "ttfc_normalization_bias_diagnostic_batch_v0.1":
            raise ValueError(f"unexpected batch schema in {path}")
        if payload.get("candidate_selected") is not None or payload.get("formal_pass_fail") is not None:
            raise ValueError("normalization-bias batch attempted selection or qualification")
        cell_id = payload["cell_id"]
        if cell_id not in rows_by_cell:
            raise ValueError(f"unknown cell {cell_id}")
        for row in payload["rows"]:
            rep = int(row["replicate"])
            if rep in rows_by_cell[cell_id]:
                raise ValueError(f"duplicate replicate {rep} in {cell_id}")
            rows_by_cell[cell_id][rep] = row

    expected = set(range(expected_n))
    out = {}
    for cell_id, rep_map in rows_by_cell.items():
        if set(rep_map) != expected:
            missing = sorted(expected - set(rep_map))
            extra = sorted(set(rep_map) - expected)
            raise ValueError(
                f"replicate coverage mismatch for {cell_id}: missing={missing[:10]} extra={extra[:10]}"
            )
        rows = [rep_map[i] for i in range(expected_n)]
        out[cell_id] = {
            key: np.asarray([float(row[key]) for row in rows], dtype=float)
            for key in WORLD_KEYS
        }
    return out


def _mechanism_for_cell(summary: dict) -> dict:
    norm = summary["endpoints"]["conditional_mean_normalized_C"]
    unscaled = summary["endpoints"]["conditional_mean_unscaled_C"]
    difference = summary["endpoints"]["conditional_mean_normalized_minus_unscaled_C"]
    normalized_bias = bool(norm["low95"] > 0.0)
    contribution = bool(difference["low95"] > 0.0)
    recentered = bool(unscaled["low95"] <= 0.0 <= unscaled["high95"])
    return {
        "normalized_bias_replication": {
            "supported": normalized_bias,
            "conditional_mean_normalized_C": norm,
        },
        "normalization_contribution": {
            "supported": contribution,
            "normalized_minus_unscaled_C": difference,
        },
        "unscaled_recentered": {
            "supported": recentered,
            "conditional_mean_unscaled_C": unscaled,
        },
        "joint_supported": bool(normalized_bias and contribution and recentered),
    }


def aggregate(rule: dict, batch_paths: list[Path]) -> dict:
    if rule.get("status") != "frozen_before_normalization_bias_outcomes":
        raise ValueError("normalization-bias rule is not frozen")
    if rule["terminal_rule"].get("candidate_selected") is not None:
        raise ValueError("diagnostic must forbid candidate selection")
    if rule["terminal_rule"].get("formal_pass_fail") is not None:
        raise ValueError("diagnostic must forbid formal pass/fail")

    cells = _collect(rule, batch_paths)
    cfg = rule["fresh_worlds"]
    n_boot = int(cfg["aggregate_bootstrap_resamples"])
    seed = int(cfg["aggregate_bootstrap_seed"])

    summaries = {}
    for i, cell in enumerate(cfg["cells"]):
        cid = cell["id"]
        summaries[cid] = {
            "geometry_mode": cell["geometry_mode"],
            "role": cell["role"],
            "n_geometry_worlds": int(len(cells[cid]["conditional_mean_normalized_C"])),
            "endpoints": {
                key: _bootstrap_mean(
                    cells[cid][key],
                    rng=np.random.default_rng(seed + 1000 * i + j),
                    n_bootstrap=n_boot,
                )
                for j, key in enumerate(WORLD_KEYS)
            },
        }

    matched = _mechanism_for_cell(summaries["DNB_MATCHED"])
    shifted = _mechanism_for_cell(summaries["DNB_SHIFTED"])

    return {
        "schema": "ttfc_normalization_bias_diagnostic_result_v0.1",
        "status": "diagnostic_complete_nonqualifying",
        "candidate_selected": None,
        "formal_pass_fail": None,
        "cell_summaries": summaries,
        "mechanism_tests": {
            "matched_normalized_bias_replication": matched["normalized_bias_replication"],
            "matched_normalization_contribution": matched["normalization_contribution"],
            "matched_unscaled_recentered": matched["unscaled_recentered"],
            "primary_normalization_bias_mechanism": {
                "supported": matched["joint_supported"],
                "requires": [
                    "matched_normalized_bias_replication",
                    "matched_normalization_contribution",
                    "matched_unscaled_recentered",
                ],
            },
            "shifted_normalized_bias_replication": shifted["normalized_bias_replication"],
            "shifted_normalization_contribution": shifted["normalization_contribution"],
            "shifted_unscaled_recentered": shifted["unscaled_recentered"],
            "shifted_joint_normalization_bias_mechanism": {
                "supported": shifted["joint_supported"],
                "requires": [
                    "shifted_normalized_bias_replication",
                    "shifted_normalization_contribution",
                    "shifted_unscaled_recentered",
                ],
            },
            "homogeneous_negative_control": {
                "normalized": summaries["DNB_HOMOGENEOUS"]["endpoints"]["conditional_mean_normalized_C"],
                "unscaled": summaries["DNB_HOMOGENEOUS"]["endpoints"]["conditional_mean_unscaled_C"],
                "normalized_minus_unscaled": summaries["DNB_HOMOGENEOUS"]["endpoints"]["conditional_mean_normalized_minus_unscaled_C"],
            },
        },
        "global_gate_predeclared": False,
        "interpretation": "This diagnostic isolates the final per-system residual-SD normalization as a possible contributor to the supported geometry-conditioned TTF-C null bias. It does not select the unscaled ablation or any other correction.",
        "terminal_rule": rule["terminal_rule"],
        "claim_firewall": rule["claim_firewall"],
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--rule", required=True)
    parser.add_argument("--batch-dir", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    rule = json.loads(Path(args.rule).read_text(encoding="utf-8"))
    paths = sorted(Path(args.batch_dir).rglob("*.json"))
    if not paths:
        raise ValueError("no normalization-bias batches found")
    result = aggregate(rule, paths)
    out = Path(args.output)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(
        json.dumps(
            {
                "primary_normalization_bias_mechanism": result["mechanism_tests"]["primary_normalization_bias_mechanism"]["supported"],
                "shifted_joint_normalization_bias_mechanism": result["mechanism_tests"]["shifted_joint_normalization_bias_mechanism"]["supported"],
                "candidate_selected": None,
                "formal_pass_fail": None,
            },
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()
