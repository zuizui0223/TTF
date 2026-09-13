from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

WORLD_KEYS = (
    "conditional_mean_C",
    "within_geometry_phase_sd_C",
    "conditional_mean_mc_se_C",
    "phase_fraction_positive_C",
    "phase_q05_C",
    "phase_q50_C",
    "phase_q95_C",
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


def _bootstrap_independent_difference(a, b, *, rng, n_bootstrap: int) -> dict:
    left = np.asarray(a, dtype=float)
    right = np.asarray(b, dtype=float)
    if left.ndim != 1 or right.ndim != 1 or len(left) < 2 or len(right) < 2:
        raise ValueError("difference bootstrap requires two vectors")
    if not np.isfinite(left).all() or not np.isfinite(right).all():
        raise ValueError("difference bootstrap inputs must be finite")
    draws = np.empty(int(n_bootstrap), dtype=float)
    for i in range(int(n_bootstrap)):
        ia = rng.integers(0, len(left), size=len(left))
        ib = rng.integers(0, len(right), size=len(right))
        draws[i] = float(np.mean(left[ia]) - np.mean(right[ib]))
    return _interval(draws, float(np.mean(left) - np.mean(right)))


def _collect(rule: dict, paths: list[Path]) -> dict[str, dict[str, np.ndarray]]:
    expected_n = int(rule["diagnostic_worlds"]["geometry_worlds_per_cell"])
    cell_rule = {cell["id"]: cell for cell in rule["cells"]}
    rows_by_cell = {cell_id: {} for cell_id in cell_rule}
    for path in paths:
        payload = json.loads(path.read_text(encoding="utf-8"))
        if payload.get("schema") != "ttfc_conditional_null_centering_batch_v0.1":
            raise ValueError(f"unexpected batch schema in {path}")
        if payload.get("candidate_selected") is not None or payload.get("formal_pass_fail") is not None:
            raise ValueError("conditional-null batch attempted selection or qualification")
        cell_id = payload["cell_id"]
        if cell_id not in rows_by_cell:
            raise ValueError(f"unknown cell {cell_id}")
        for row in payload["rows"]:
            rep = int(row["replicate"])
            if rep in rows_by_cell[cell_id]:
                raise ValueError(f"duplicate replicate {rep} in {cell_id}")
            rows_by_cell[cell_id][rep] = row

    out = {}
    expected = set(range(expected_n))
    for cell_id, rep_map in rows_by_cell.items():
        if set(rep_map) != expected:
            missing = sorted(expected - set(rep_map))
            extra = sorted(set(rep_map) - expected)
            raise ValueError(f"replicate coverage mismatch for {cell_id}: missing={missing[:10]} extra={extra[:10]}")
        rows = [rep_map[i] for i in range(expected_n)]
        out[cell_id] = {
            key: np.asarray([float(row[key]) for row in rows], dtype=float)
            for key in WORLD_KEYS
        }
    return out


def aggregate(rule: dict, batch_paths: list[Path]) -> dict:
    if rule.get("status") != "frozen_before_conditional_null_centering_outcomes":
        raise ValueError("conditional-null centering rule is not frozen")
    if rule["terminal_rule"].get("candidate_selected") is not None:
        raise ValueError("diagnostic must forbid candidate selection")
    if rule["terminal_rule"].get("formal_pass_fail") is not None:
        raise ValueError("diagnostic must forbid formal pass/fail")

    cells = _collect(rule, batch_paths)
    cfg = rule["diagnostic_worlds"]
    n_boot = int(cfg["aggregate_bootstrap_resamples"])
    seed = int(cfg["aggregate_bootstrap_seed"])

    summaries = {}
    for i, cell in enumerate(rule["cells"]):
        cid = cell["id"]
        summaries[cid] = {
            "geometry_mode": cell["geometry_mode"],
            "n_geometry_worlds": int(len(cells[cid]["conditional_mean_C"])),
            "endpoints": {
                key: _bootstrap_mean(
                    cells[cid][key],
                    rng=np.random.default_rng(seed + 1000 * i + j),
                    n_bootstrap=n_boot,
                )
                for j, key in enumerate(WORLD_KEYS)
            },
        }

    hom = cells["DCN_HOMOGENEOUS"]["conditional_mean_C"]
    matched = cells["DCN_MATCHED"]["conditional_mean_C"]
    shifted = cells["DCN_SHIFTED"]["conditional_mean_C"]

    matched_mean = summaries["DCN_MATCHED"]["endpoints"]["conditional_mean_C"]
    shifted_mean = summaries["DCN_SHIFTED"]["endpoints"]["conditional_mean_C"]
    matched_bias = bool(matched_mean["low95"] > 0.0)
    shifted_bias = bool(shifted_mean["low95"] > 0.0)

    matched_minus_hom = _bootstrap_independent_difference(
        matched,
        hom,
        rng=np.random.default_rng(seed + 100_000),
        n_bootstrap=n_boot,
    )
    shifted_minus_hom = _bootstrap_independent_difference(
        shifted,
        hom,
        rng=np.random.default_rng(seed + 100_001),
        n_bootstrap=n_boot,
    )
    shifted_minus_matched = _bootstrap_independent_difference(
        shifted,
        matched,
        rng=np.random.default_rng(seed + 100_002),
        n_bootstrap=n_boot,
    )
    matched_amp = bool(matched_minus_hom["low95"] > 0.0)
    shifted_amp = bool(shifted_minus_hom["low95"] > 0.0)

    return {
        "schema": "ttfc_conditional_null_centering_result_v0.1",
        "status": "diagnostic_complete_nonqualifying",
        "candidate_selected": None,
        "formal_pass_fail": None,
        "cell_summaries": summaries,
        "mechanism_tests": {
            "matched_conditional_null_bias": {
                "supported": matched_bias,
                "conditional_mean_C": matched_mean,
            },
            "shifted_conditional_null_bias": {
                "supported": shifted_bias,
                "conditional_mean_C": shifted_mean,
            },
            "matched_heterogeneity_amplification": {
                "supported": matched_amp,
                "matched_minus_homogeneous": matched_minus_hom,
            },
            "shifted_heterogeneity_amplification": {
                "supported": shifted_amp,
                "shifted_minus_homogeneous": shifted_minus_hom,
            },
            "shifted_minus_matched": shifted_minus_matched,
            "primary_geometry_conditioned_null_bias": {
                "supported": bool(matched_bias and matched_amp),
                "requires": [
                    "matched_conditional_null_bias",
                    "matched_heterogeneity_amplification",
                ],
            },
        },
        "global_gate_predeclared": False,
        "interpretation": "This diagnostic tests geometry-conditioned null centering only; it does not select a correction or qualify TTF-C.",
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
        raise ValueError("no conditional-null batches found")
    result = aggregate(rule, paths)
    out = Path(args.output)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps({
        "primary_geometry_conditioned_null_bias": result["mechanism_tests"]["primary_geometry_conditioned_null_bias"]["supported"],
        "matched_conditional_null_bias": result["mechanism_tests"]["matched_conditional_null_bias"]["supported"],
        "matched_heterogeneity_amplification": result["mechanism_tests"]["matched_heterogeneity_amplification"]["supported"],
        "candidate_selected": None,
        "formal_pass_fail": None,
    }, sort_keys=True))


if __name__ == "__main__":
    main()
