from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

from ttf.core import spearman_rho


WORLD_KEYS = (
    "raw_M",
    "raw_C",
    "loo_score_sd_M",
    "loo_score_sd_C",
    "loo_score_mean_M",
    "loo_score_mean_C",
    "loo_score_min_M",
    "loo_score_min_C",
    "mean_abs_loo_score_shift_M",
    "mean_abs_loo_score_shift_C",
    "mean_prediction_loo_sd_M",
    "mean_prediction_loo_sd_C",
    "mean_effective_training_system_count",
    "mean_dominant_training_system_share",
)


def _interval(samples: np.ndarray, estimate: float) -> dict:
    x = np.asarray(samples, dtype=float)
    if x.ndim != 1 or len(x) < 1 or not np.isfinite(x).all():
        raise ValueError("bootstrap samples must be finite 1D")
    return {
        "estimate": float(estimate),
        "low95": float(np.quantile(x, 0.025)),
        "high95": float(np.quantile(x, 0.975)),
    }


def _bootstrap_mean(values: np.ndarray, *, rng: np.random.Generator, n_bootstrap: int) -> dict:
    x = np.asarray(values, dtype=float)
    if x.ndim != 1 or len(x) < 2 or not np.isfinite(x).all():
        raise ValueError("mean bootstrap requires finite 1D values")
    draws = np.empty(int(n_bootstrap), dtype=float)
    for i in range(int(n_bootstrap)):
        idx = rng.integers(0, len(x), size=len(x))
        draws[i] = float(x[idx].mean())
    return _interval(draws, float(x.mean()))


def _bootstrap_spearman(x: np.ndarray, y: np.ndarray, *, rng: np.random.Generator, n_bootstrap: int) -> dict:
    xx = np.asarray(x, dtype=float)
    yy = np.asarray(y, dtype=float)
    if xx.ndim != 1 or yy.shape != xx.shape or len(xx) < 3:
        raise ValueError("Spearman bootstrap requires equal 1D vectors")
    if not np.isfinite(xx).all() or not np.isfinite(yy).all():
        raise ValueError("Spearman bootstrap inputs must be finite")
    estimate = float(spearman_rho(xx, yy))
    draws = np.empty(int(n_bootstrap), dtype=float)
    for i in range(int(n_bootstrap)):
        idx = rng.integers(0, len(xx), size=len(xx))
        draws[i] = float(spearman_rho(xx[idx], yy[idx]))
    return _interval(draws, estimate)


def _bootstrap_independent_mean_difference(left: np.ndarray, right: np.ndarray, *, rng: np.random.Generator, n_bootstrap: int) -> dict:
    a = np.asarray(left, dtype=float)
    b = np.asarray(right, dtype=float)
    if a.ndim != 1 or b.ndim != 1 or len(a) < 2 or len(b) < 2:
        raise ValueError("independent difference needs two vectors")
    if not np.isfinite(a).all() or not np.isfinite(b).all():
        raise ValueError("independent difference inputs must be finite")
    draws = np.empty(int(n_bootstrap), dtype=float)
    for i in range(int(n_bootstrap)):
        ia = rng.integers(0, len(a), size=len(a))
        ib = rng.integers(0, len(b), size=len(b))
        draws[i] = float(a[ia].mean() - b[ib].mean())
    return _interval(draws, float(a.mean() - b.mean()))


def _collect(rule: dict, batch_paths: list[Path]) -> dict[str, dict[str, np.ndarray]]:
    expected_n = int(rule["diagnostic_worlds"]["worlds_per_cell"])
    cell_rule = {cell["id"]: cell for cell in rule["cells"]}
    rows_by_cell: dict[str, dict[int, dict]] = {cell_id: {} for cell_id in cell_rule}

    for path in batch_paths:
        payload = json.loads(path.read_text(encoding="utf-8"))
        if payload.get("schema") != "ttfc_training_leverage_diagnostic_batch_v0.1":
            raise ValueError(f"unexpected batch schema in {path}")
        if payload.get("candidate_selected") is not None or payload.get("formal_pass_fail") is not None:
            raise ValueError("training-leverage batch attempted a forbidden decision")
        cell_id = payload["cell_id"]
        if cell_id not in rows_by_cell:
            raise ValueError(f"unknown cell {cell_id}")
        for row in payload["rows"]:
            rep = int(row["replicate"])
            if rep in rows_by_cell[cell_id]:
                raise ValueError(f"duplicate replicate {rep} in {cell_id}")
            if set(row["world"]) != set(WORLD_KEYS):
                raise ValueError(f"world endpoint drift in {cell_id} replicate {rep}")
            rows_by_cell[cell_id][rep] = row

    out: dict[str, dict[str, np.ndarray]] = {}
    expected = set(range(expected_n))
    for cell_id, rep_map in rows_by_cell.items():
        if set(rep_map) != expected:
            missing = sorted(expected - set(rep_map))
            extra = sorted(set(rep_map) - expected)
            raise ValueError(f"replicate coverage mismatch for {cell_id}: missing={missing[:10]} extra={extra[:10]}")
        rows = [rep_map[i] for i in range(expected_n)]
        out[cell_id] = {
            key: np.asarray([float(row["world"][key]) for row in rows], dtype=float)
            for key in WORLD_KEYS
        }
    return out


def _field_leverage_test(data: dict[str, np.ndarray], *, base_seed: int, offset: int, n_bootstrap: int) -> dict:
    score_sd = _bootstrap_spearman(
        data["raw_C"], data["loo_score_sd_C"],
        rng=np.random.default_rng(base_seed + offset), n_bootstrap=n_bootstrap,
    )
    pred_sd = _bootstrap_spearman(
        data["raw_C"], data["mean_prediction_loo_sd_C"],
        rng=np.random.default_rng(base_seed + offset + 1), n_bootstrap=n_bootstrap,
    )
    return {
        "raw_C_vs_loo_score_sd_C_spearman": score_sd,
        "raw_C_vs_mean_prediction_loo_sd_C_spearman": pred_sd,
        "supported": bool(score_sd["low95"] > 0.0 and pred_sd["low95"] > 0.0),
    }


def aggregate(rule: dict, batch_paths: list[Path]) -> dict:
    if rule.get("status") != "frozen_before_training_leverage_diagnostic_outcomes":
        raise ValueError("training-leverage diagnostic rule is not frozen")
    if rule["terminal_rule"].get("candidate_selected") is not None:
        raise ValueError("diagnostic rule must forbid candidate selection")
    if rule["terminal_rule"].get("formal_pass_fail") is not None:
        raise ValueError("diagnostic rule must forbid formal pass/fail")

    cells = _collect(rule, batch_paths)
    cfg = rule["diagnostic_worlds"]
    n_bootstrap = int(cfg["aggregate_bootstrap_resamples"])
    base_seed = int(cfg["aggregate_bootstrap_seed"])

    cell_summaries = {}
    for i, cell in enumerate(rule["cells"]):
        data = cells[cell["id"]]
        endpoints = {}
        for j, key in enumerate(WORLD_KEYS):
            endpoints[key] = _bootstrap_mean(
                data[key],
                rng=np.random.default_rng(base_seed + 1000 * i + j),
                n_bootstrap=n_bootstrap,
            )
        cell_summaries[cell["id"]] = {
            "response_mode": cell["response_mode"],
            "geometry_profile": cell["geometry_profile"],
            "n_worlds": int(len(data["raw_C"])),
            "endpoints": endpoints,
        }

    private_c = cells["DTL_PRIVATE_C"]
    private_m = cells["DTL_PRIVATE_M"]
    private_c_field = _field_leverage_test(private_c, base_seed=base_seed, offset=100_000, n_bootstrap=n_bootstrap)
    private_m_field = _field_leverage_test(private_m, base_seed=base_seed, offset=100_100, n_bootstrap=n_bootstrap)

    dom = _bootstrap_spearman(
        private_c["raw_C"], private_c["mean_dominant_training_system_share"],
        rng=np.random.default_rng(base_seed + 110_000), n_bootstrap=n_bootstrap,
    )
    eff = _bootstrap_spearman(
        private_c["raw_C"], private_c["mean_effective_training_system_count"],
        rng=np.random.default_rng(base_seed + 110_001), n_bootstrap=n_bootstrap,
    )
    private_c_concentration = {
        "raw_C_vs_dominant_share_spearman": dom,
        "raw_C_vs_effective_system_count_spearman": eff,
        "supported": bool(dom["low95"] > 0.0 and eff["high95"] < 0.0),
    }

    matched = cells["DTL_POWER_M_MATCHED"]
    shifted = cells["DTL_POWER_M_SHIFTED"]
    shift_instability = _bootstrap_independent_mean_difference(
        shifted["loo_score_sd_C"], matched["loo_score_sd_C"],
        rng=np.random.default_rng(base_seed + 120_000), n_bootstrap=n_bootstrap,
    )
    shift_raw = _bootstrap_independent_mean_difference(
        shifted["raw_C"], matched["raw_C"],
        rng=np.random.default_rng(base_seed + 120_001), n_bootstrap=n_bootstrap,
    )
    shifted_test = {
        "loo_score_sd_C_shifted_minus_matched": shift_instability,
        "raw_C_shifted_minus_matched": shift_raw,
        "comparison": "independent bootstrap because matched and shifted cells use independently generated worlds",
        "supported": bool(shift_instability["low95"] > 0.0 and shift_raw["high95"] < 0.0),
    }

    specificity = {}
    specificity_supported = False
    for j, cell_id in enumerate(("DTL_PRIVATE_M", "DTL_POWER_M_MATCHED")):
        data = cells[cell_id]
        interval = _bootstrap_mean(
            data["loo_score_sd_C"] - data["loo_score_sd_M"],
            rng=np.random.default_rng(base_seed + 130_000 + j),
            n_bootstrap=n_bootstrap,
        )
        flag = bool(interval["low95"] > 0.0)
        specificity[cell_id] = {"C_minus_M_loo_score_sd": interval, "supported": flag}
        specificity_supported = specificity_supported or flag
    specificity["supported"] = bool(specificity_supported)

    negative_control = {}
    for cell_id in ("DTL_POWER_C_MATCHED", "DTL_POWER_C_SHIFTED"):
        e = cell_summaries[cell_id]["endpoints"]
        negative_control[cell_id] = {
            "raw_C": e["raw_C"],
            "loo_score_mean_C": e["loo_score_mean_C"],
            "loo_score_min_C": e["loo_score_min_C"],
            "loo_score_sd_C": e["loo_score_sd_C"],
        }

    return {
        "schema": "ttfc_training_leverage_diagnostic_result_v0.1",
        "status": "diagnostic_complete_nonqualifying",
        "candidate_selected": None,
        "formal_pass_fail": None,
        "cell_summaries": cell_summaries,
        "mechanism_tests": {
            "private_relation_field_leverage": private_c_field,
            "private_relation_geometry_concentration": private_c_concentration,
            "private_mismatch_field_leverage": private_m_field,
            "shifted_power_training_instability": shifted_test,
            "estimand_specificity": specificity,
            "shared_relation_negative_control": negative_control,
        },
        "global_gate_predeclared": False,
        "interpretation": "Report each predeclared training-leverage mechanism test separately; no repair candidate or TTF-C qualification is selected.",
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
    batch_paths = sorted(Path(args.batch_dir).rglob("*.json"))
    if not batch_paths:
        raise ValueError("no training-leverage diagnostic batches found")
    result = aggregate(rule, batch_paths)
    out = Path(args.output)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps({
        "private_relation_field_leverage": result["mechanism_tests"]["private_relation_field_leverage"]["supported"],
        "private_relation_geometry_concentration": result["mechanism_tests"]["private_relation_geometry_concentration"]["supported"],
        "shifted_power_training_instability": result["mechanism_tests"]["shifted_power_training_instability"]["supported"],
        "estimand_specificity": result["mechanism_tests"]["estimand_specificity"]["supported"],
        "candidate_selected": None,
        "formal_pass_fail": None,
    }, sort_keys=True))


if __name__ == "__main__":
    main()
