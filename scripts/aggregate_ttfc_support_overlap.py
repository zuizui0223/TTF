from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

from ttf.core import spearman_rho


WORLD_KEYS = (
    "raw_M",
    "partial_M",
    "alignment_M",
    "raw_C",
    "partial_C",
    "alignment_C",
    "mean_opportunity",
    "mean_prior_fraction",
    "low_support_fraction",
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


def _bootstrap_mean(
    values: np.ndarray,
    *,
    rng: np.random.Generator,
    n_bootstrap: int,
) -> dict:
    x = np.asarray(values, dtype=float)
    if x.ndim != 1 or len(x) < 2 or not np.isfinite(x).all():
        raise ValueError("mean bootstrap requires finite 1D values")
    draws = np.empty(int(n_bootstrap), dtype=float)
    for i in range(int(n_bootstrap)):
        idx = rng.integers(0, len(x), size=len(x))
        draws[i] = float(x[idx].mean())
    return _interval(draws, float(x.mean()))


def _bootstrap_spearman(
    x: np.ndarray,
    y: np.ndarray,
    *,
    rng: np.random.Generator,
    n_bootstrap: int,
) -> dict:
    xx = np.asarray(x, dtype=float)
    yy = np.asarray(y, dtype=float)
    if xx.ndim != 1 or yy.shape != xx.shape or len(xx) < 3:
        raise ValueError("Spearman bootstrap requires equal finite vectors")
    if not np.isfinite(xx).all() or not np.isfinite(yy).all():
        raise ValueError("Spearman bootstrap inputs must be finite")
    estimate = float(spearman_rho(xx, yy))
    draws = np.empty(int(n_bootstrap), dtype=float)
    for i in range(int(n_bootstrap)):
        idx = rng.integers(0, len(xx), size=len(xx))
        draws[i] = float(spearman_rho(xx[idx], yy[idx]))
    return _interval(draws, estimate)


def _bootstrap_independent_mean_difference(
    left: np.ndarray,
    right: np.ndarray,
    *,
    rng: np.random.Generator,
    n_bootstrap: int,
) -> dict:
    """Bootstrap mean(left)-mean(right) for independently generated cells."""
    a = np.asarray(left, dtype=float)
    b = np.asarray(right, dtype=float)
    if a.ndim != 1 or b.ndim != 1 or len(a) < 2 or len(b) < 2:
        raise ValueError("independent difference needs two finite vectors")
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
        if payload.get("schema") != "ttfc_support_overlap_diagnostic_batch_v0.1":
            raise ValueError(f"unexpected batch schema in {path}")
        if payload.get("candidate_selected") is not None:
            raise ValueError("support diagnostic batch attempted candidate selection")
        if payload.get("formal_pass_fail") is not None:
            raise ValueError("support diagnostic batch attempted formal decision")
        cell_id = payload["cell_id"]
        if cell_id not in rows_by_cell:
            raise ValueError(f"unknown diagnostic cell: {cell_id}")
        for row in payload["rows"]:
            rep = int(row["replicate"])
            if rep in rows_by_cell[cell_id]:
                raise ValueError(f"duplicate replicate {rep} in {cell_id}")
            world = row["world"]
            if set(world) != set(WORLD_KEYS):
                raise ValueError(f"world endpoint drift in {cell_id} replicate {rep}")
            rows_by_cell[cell_id][rep] = row

    out: dict[str, dict[str, np.ndarray]] = {}
    expected = set(range(expected_n))
    for cell_id, rep_map in rows_by_cell.items():
        if set(rep_map) != expected:
            missing = sorted(expected - set(rep_map))
            extra = sorted(set(rep_map) - expected)
            raise ValueError(
                f"replicate coverage mismatch for {cell_id}: "
                f"missing={missing[:10]} extra={extra[:10]}"
            )
        rows = [rep_map[i] for i in range(expected_n)]
        out[cell_id] = {
            key: np.asarray([float(row["world"][key]) for row in rows], dtype=float)
            for key in WORLD_KEYS
        }
    return out


def aggregate(rule: dict, batch_paths: list[Path]) -> dict:
    if rule.get("status") != "frozen_before_support_overlap_diagnostic_outcomes":
        raise ValueError("support-overlap diagnostic rule is not frozen")
    if rule["terminal_rule"].get("candidate_selected") is not None:
        raise ValueError("diagnostic terminal rule must forbid candidate selection")
    if rule["terminal_rule"].get("formal_pass_fail") is not None:
        raise ValueError("diagnostic terminal rule must forbid formal pass/fail")

    cells = _collect(rule, batch_paths)
    cfg = rule["diagnostic_worlds"]
    n_bootstrap = int(cfg["aggregate_bootstrap_resamples"])
    base_seed = int(cfg["aggregate_bootstrap_seed"])

    cell_summaries = {}
    for index, cell in enumerate(rule["cells"]):
        cell_id = cell["id"]
        data = cells[cell_id]
        endpoint_summary = {}
        for j, key in enumerate(WORLD_KEYS):
            rng = np.random.default_rng(base_seed + 1000 * index + j)
            endpoint_summary[key] = _bootstrap_mean(
                data[key], rng=rng, n_bootstrap=n_bootstrap
            )
        cell_summaries[cell_id] = {
            "response_mode": cell["response_mode"],
            "geometry_profile": cell["geometry_profile"],
            "n_worlds": int(len(data["raw_C"])),
            "endpoints": endpoint_summary,
        }

    private_tests = {}
    for offset, cell_id in enumerate(("DSO_PRIVATE_M", "DSO_PRIVATE_C")):
        data = cells[cell_id]
        delta_c = data["raw_C"] - data["partial_C"]
        delta_m = data["raw_M"] - data["partial_M"]
        private_tests[cell_id] = {
            "coupling_raw_minus_partial": _bootstrap_mean(
                delta_c,
                rng=np.random.default_rng(base_seed + 100_000 + 10 * offset),
                n_bootstrap=n_bootstrap,
            ),
            "mismatch_raw_minus_partial": _bootstrap_mean(
                delta_m,
                rng=np.random.default_rng(base_seed + 100_001 + 10 * offset),
                n_bootstrap=n_bootstrap,
            ),
            "raw_C_vs_alignment_C_spearman": _bootstrap_spearman(
                data["raw_C"],
                data["alignment_C"],
                rng=np.random.default_rng(base_seed + 100_002 + 10 * offset),
                n_bootstrap=n_bootstrap,
            ),
            "raw_M_vs_alignment_M_spearman": _bootstrap_spearman(
                data["raw_M"],
                data["alignment_M"],
                rng=np.random.default_rng(base_seed + 100_003 + 10 * offset),
                n_bootstrap=n_bootstrap,
            ),
            "C_minus_M_mediation_effect": _bootstrap_mean(
                delta_c - delta_m,
                rng=np.random.default_rng(base_seed + 100_004 + 10 * offset),
                n_bootstrap=n_bootstrap,
            ),
        }

    private_relation = private_tests["DSO_PRIVATE_C"]
    private_relation_supported = bool(
        private_relation["coupling_raw_minus_partial"]["low95"] > 0.0
        and private_relation["raw_C_vs_alignment_C_spearman"]["low95"] > 0.0
    )
    private_mismatch = private_tests["DSO_PRIVATE_M"]
    private_mismatch_supported = bool(
        private_mismatch["coupling_raw_minus_partial"]["low95"] > 0.0
        and private_mismatch["raw_C_vs_alignment_C_spearman"]["low95"] > 0.0
    )

    matched = cells["DSO_POWER_M_MATCHED"]
    shifted = cells["DSO_POWER_M_SHIFTED"]
    prior_shift = _bootstrap_independent_mean_difference(
        shifted["mean_prior_fraction"],
        matched["mean_prior_fraction"],
        rng=np.random.default_rng(base_seed + 200_000),
        n_bootstrap=n_bootstrap,
    )
    raw_c_shift = _bootstrap_independent_mean_difference(
        shifted["raw_C"],
        matched["raw_C"],
        rng=np.random.default_rng(base_seed + 200_001),
        n_bootstrap=n_bootstrap,
    )
    shifted_power_support_loss_supported = bool(
        prior_shift["low95"] > 0.0 and raw_c_shift["high95"] < 0.0
    )

    specificity_private_m = bool(
        private_tests["DSO_PRIVATE_M"]["C_minus_M_mediation_effect"]["low95"] > 0.0
    )
    specificity_private_c = bool(
        private_tests["DSO_PRIVATE_C"]["C_minus_M_mediation_effect"]["low95"] > 0.0
    )
    estimand_specificity_supported = bool(
        specificity_private_m or specificity_private_c
    )

    negative_control = {}
    for cell_id in ("DSO_POWER_C_MATCHED", "DSO_POWER_C_SHIFTED"):
        data = cells[cell_id]
        negative_control[cell_id] = {
            "raw_C": cell_summaries[cell_id]["endpoints"]["raw_C"],
            "partial_C": cell_summaries[cell_id]["endpoints"]["partial_C"],
            "raw_minus_partial_C": _bootstrap_mean(
                data["raw_C"] - data["partial_C"],
                rng=np.random.default_rng(
                    base_seed + 300_000 + (0 if cell_id.endswith("MATCHED") else 1)
                ),
                n_bootstrap=n_bootstrap,
            ),
        }

    return {
        "schema": "ttfc_support_overlap_diagnostic_result_v0.1",
        "status": "diagnostic_complete_nonqualifying",
        "candidate_selected": None,
        "formal_pass_fail": None,
        "cell_summaries": cell_summaries,
        "mechanism_tests": {
            "private_relation_opportunity_mediation": {
                "supported": private_relation_supported,
                **private_relation,
            },
            "private_mismatch_opportunity_mediation": {
                "supported": private_mismatch_supported,
                **private_mismatch,
            },
            "shifted_power_support_loss": {
                "supported": shifted_power_support_loss_supported,
                "mean_prior_fraction_shifted_minus_matched": prior_shift,
                "raw_C_shifted_minus_matched": raw_c_shift,
                "comparison": "independent bootstrap because matched and shifted cells use independently generated worlds",
            },
            "estimand_specificity": {
                "supported": estimand_specificity_supported,
                "private_M_specificity_supported": specificity_private_m,
                "private_C_specificity_supported": specificity_private_c,
                "private_M_C_minus_M_mediation_effect": private_tests[
                    "DSO_PRIVATE_M"
                ]["C_minus_M_mediation_effect"],
                "private_C_C_minus_M_mediation_effect": private_tests[
                    "DSO_PRIVATE_C"
                ]["C_minus_M_mediation_effect"],
            },
            "shared_relation_negative_control": negative_control,
        },
        "global_gate_predeclared": False,
        "interpretation": "Report each predeclared mechanism test separately; this diagnostic does not select a repair candidate or qualify TTF-C.",
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
        raise ValueError("no support-overlap diagnostic batches found")
    result = aggregate(rule, batch_paths)
    out = Path(args.output)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(
        json.dumps(
            {
                "private_relation_supported": result["mechanism_tests"][
                    "private_relation_opportunity_mediation"
                ]["supported"],
                "shifted_power_support_loss_supported": result["mechanism_tests"][
                    "shifted_power_support_loss"
                ]["supported"],
                "estimand_specificity_supported": result["mechanism_tests"][
                    "estimand_specificity"
                ]["supported"],
                "candidate_selected": None,
                "formal_pass_fail": None,
            },
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()
