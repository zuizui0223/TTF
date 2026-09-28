#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

from ttf.relational_benchmark import known_truth_design
from ttf.relational_qualification import (
    dyadic_signal_support,
    residualize_primary_relation,
)


def one_scenario(*, f, b, g, concentration, seed, n_source, n_target):
    design = known_truth_design(
        n_source=n_source,
        n_target=n_target,
        endpoint_retained_fraction=f,
        baseline_survival_after_fe=b,
        geography_survival_after_baseline=g,
        concentration=concentration,
        seed=seed,
    )
    c1_residual, c1 = residualize_primary_relation(
        design.source,
        design.target,
        design.primary,
        design.baseline_control[:, None],
    )
    c2_residual, c2 = residualize_primary_relation(
        design.source,
        design.target,
        design.primary,
        np.column_stack((
            design.baseline_control,
            design.geographic_control,
        )),
    )
    support = dyadic_signal_support(
        design.source,
        design.target,
        c2_residual,
    )
    observed_geo = float(
        np.dot(c2_residual, c2_residual)
        / np.dot(c1_residual, c1_residual)
    )
    truth = design.truth
    return {
        "f": f,
        "b": b,
        "g": g,
        "concentration": concentration,
        "seed": seed,
        "truth": truth,
        "observed": {
            "endpoint_retained_fraction": c1.source_target_fe_retained_variance_fraction,
            "baseline_survival_after_fe": c1.control_unique_variance_fraction_after_fe,
            "geography_survival_after_baseline": observed_geo,
            "C1_total_unique_variance_fraction": c1.total_unique_variance_fraction,
            "C2_total_unique_variance_fraction": c2.total_unique_variance_fraction,
            "signal_effective_sources": support.signal_effective_sources,
            "signal_effective_targets": support.signal_effective_targets,
            "max_source_signal_share": support.max_source_signal_share,
            "max_target_signal_share": support.max_target_signal_share,
        },
        "absolute_error": {
            "endpoint_retained_fraction": abs(
                c1.source_target_fe_retained_variance_fraction
                - float(truth["endpoint_retained_fraction"])
            ),
            "baseline_survival_after_fe": abs(
                c1.control_unique_variance_fraction_after_fe
                - float(truth["baseline_survival_after_fe"])
            ),
            "geography_survival_after_baseline": abs(
                observed_geo
                - float(truth["geography_survival_after_baseline"])
            ),
            "C1_total_unique_variance_fraction": abs(
                c1.total_unique_variance_fraction
                - float(truth["C1_total_unique_variance_fraction"])
            ),
            "C2_total_unique_variance_fraction": abs(
                c2.total_unique_variance_fraction
                - float(truth["C2_total_unique_variance_fraction"])
            ),
        },
    }


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument(
        "--contract",
        type=Path,
        default=Path("docs/supporting/ttf_q_known_truth_benchmark_v0.1.json"),
    )
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()

    contract = json.loads(args.contract.read_text())
    if contract.get("schema") != "ttf_q_known_truth_benchmark_v0.1":
        raise RuntimeError("unexpected known-truth benchmark contract")
    if contract.get("status") != "FROZEN_BEFORE_BENCHMARK_RESULTS":
        raise RuntimeError("known-truth benchmark contract is not pre-result frozen")

    grid = contract["decomposition_grid"]
    rows = []
    for f in map(float, grid["endpoint_retained_fraction"]):
        for b in map(float, grid["baseline_survival_after_fe"]):
            for g in map(float, grid["geography_survival_after_baseline"]):
                for seed in map(int, grid["seeds"]):
                    for concentration in map(str, grid["concentration"]):
                        rows.append(one_scenario(
                            f=f,
                            b=b,
                            g=g,
                            concentration=concentration,
                            seed=seed,
                            n_source=int(grid["n_source"]),
                            n_target=int(grid["n_target"]),
                        ))

    fields = list(rows[0]["absolute_error"])
    max_error = {
        field: float(max(row["absolute_error"][field] for row in rows))
        for field in fields
    }

    keyed = {
        (
            row["f"], row["b"], row["g"], row["seed"], row["concentration"]
        ): row
        for row in rows
    }
    comparisons = []
    for f in map(float, grid["endpoint_retained_fraction"]):
        for b in map(float, grid["baseline_survival_after_fe"]):
            for g in map(float, grid["geography_survival_after_baseline"]):
                for seed in map(int, grid["seeds"]):
                    broad = keyed[(f, b, g, seed, "broad")]
                    source = keyed[(f, b, g, seed, "source")]
                    target = keyed[(f, b, g, seed, "target")]
                    comparisons.append({
                        "f": f, "b": b, "g": g, "seed": seed,
                        "source_concentration_detected": (
                            source["observed"]["signal_effective_sources"]
                            < broad["observed"]["signal_effective_sources"]
                        ),
                        "target_concentration_detected": (
                            target["observed"]["signal_effective_targets"]
                            < broad["observed"]["signal_effective_targets"]
                        ),
                        "source_effective_ratio": (
                            source["observed"]["signal_effective_sources"]
                            / broad["observed"]["signal_effective_sources"]
                        ),
                        "target_effective_ratio": (
                            target["observed"]["signal_effective_targets"]
                            / broad["observed"]["signal_effective_targets"]
                        ),
                    })

    payload = {
        "schema": "ttf_q_known_truth_decomposition_result_v0.1",
        "status": "COMPLETE_KNOWN_TRUTH_INFORMATION_DECOMPOSITION",
        "contract_status": contract["status"],
        "scenarios": len(rows),
        "nested_identity_max_absolute_error": max_error,
        "nested_identity_max_error_overall": float(max(max_error.values())),
        "concentration_comparisons": len(comparisons),
        "source_concentration_detection_fraction": float(np.mean([
            row["source_concentration_detected"] for row in comparisons
        ])),
        "target_concentration_detection_fraction": float(np.mean([
            row["target_concentration_detected"] for row in comparisons
        ])),
        "median_source_effective_ratio_when_concentrated": float(np.median([
            row["source_effective_ratio"] for row in comparisons
        ])),
        "median_target_effective_ratio_when_concentrated": float(np.median([
            row["target_effective_ratio"] for row in comparisons
        ])),
        "rows": rows,
        "boundary": contract["boundary"],
    }
    if payload["nested_identity_max_error_overall"] > 1e-10:
        raise RuntimeError("known-truth decomposition failed exact recovery tolerance")

    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    print(json.dumps({
        key: payload[key]
        for key in (
            "status",
            "scenarios",
            "nested_identity_max_error_overall",
            "source_concentration_detection_fraction",
            "target_concentration_detection_fraction",
            "median_source_effective_ratio_when_concentrated",
            "median_target_effective_ratio_when_concentrated",
        )
    }, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
