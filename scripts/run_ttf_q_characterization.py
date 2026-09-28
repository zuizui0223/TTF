#!/usr/bin/env python3
from __future__ import annotations

import argparse
from dataclasses import asdict
import hashlib
import json
from pathlib import Path

import numpy as np

from ttf.relational_qualification import (
    calibrated_detectability_envelope,
    detectability_surface,
    dyadic_signal_support,
    minimum_detectable_effects,
    relation_repeatability,
    residualize_primary_relation,
)


def sha256_path(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def zscore(values: np.ndarray) -> np.ndarray:
    x = np.asarray(values, dtype=float)
    sd = float(x.std())
    if not np.isfinite(x).all() or sd <= np.finfo(float).eps:
        raise ValueError("non-finite or constant predictor")
    return (x - float(x.mean())) / sd


def remap(values: np.ndarray) -> np.ndarray:
    unique = sorted(map(int, np.unique(values)))
    lookup = {value: i for i, value in enumerate(unique)}
    return np.asarray([lookup[int(value)] for value in values], dtype=np.int64)


def build_design(data: np.lib.npyio.NpzFile, kind: str, panel: str):
    source = np.asarray(data[f"{panel}_source_index"], np.int64)
    target = np.asarray(data[f"{panel}_target_index"], np.int64)
    coverage = np.asarray(data[f"{panel}_coverage"], float)
    same_class = np.asarray(data[f"{panel}_same_class"], float)
    same_order = np.asarray(data[f"{panel}_same_order"], float)
    same_family = np.asarray(data[f"{panel}_same_family"], float)
    locality_ratio = np.asarray(data[f"{panel}_locality_count_ratio"], float)

    common = [
        zscore(coverage),
        same_class - same_class.mean(),
        same_order - same_order.mean(),
        same_family - same_family.mean(),
        zscore(np.abs(np.log(locality_ratio))),
    ]
    if kind in {"environment", "current"}:
        relation_key = (
            f"{panel}_R_env"
            if kind == "environment"
            else f"{panel}_R_current"
        )
        primary = zscore(np.asarray(data[relation_key], float))
        names = [
            "z_R_env" if kind == "environment" else "z_R_current",
            "z_geographic_coverage",
            "centered_same_class",
            "centered_same_order",
            "centered_same_family",
            "z_abs_log_locality_count_ratio",
        ]
        predictors = np.column_stack([primary, *common])
    else:
        primary = zscore(np.asarray(data[f"{panel}_R_hist"], float))
        current = zscore(np.asarray(data[f"{panel}_R_current"], float))
        names = [
            "z_R_hist",
            "z_R_current",
            "z_geographic_coverage",
            "centered_same_class",
            "centered_same_order",
            "centered_same_family",
            "z_abs_log_locality_count_ratio",
        ]
        predictors = np.column_stack([primary, current, *common])
    return remap(source), remap(target), names, predictors


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--design", type=Path, required=True)
    parser.add_argument(
        "--contract",
        type=Path,
        default=Path("docs/supporting/ttf_q_v0.1.json"),
    )
    parser.add_argument(
        "--kind",
        choices=("environment", "current", "historical"),
        required=True,
    )
    parser.add_argument(
        "--panel",
        choices=("development", "confirmatory"),
        default="development",
    )
    parser.add_argument("--relation-replicates", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--worlds-per-cell", type=int)
    args = parser.parse_args()

    contract = json.loads(args.contract.read_text())
    if contract.get("schema") != "ttf_q_response_blind_qualification_v0.1":
        raise RuntimeError("unexpected TTF-Q contract")
    if contract["boundary"]["closed_v03_final_state"] != "CLOSED_NO_EVALUABLE_TEST":
        raise RuntimeError("TTF-Q cannot rewrite the frozen v0.3 terminal state")
    if contract["boundary"]["genetic_response_opening_authorized"]:
        raise RuntimeError("TTF-Q must remain response blind")

    data = np.load(args.design, allow_pickle=False)
    source, target, predictor_names, predictors = build_design(
        data, args.kind, args.panel
    )
    primary = predictors[:, 0]
    controls = predictors[:, 1:]

    residual, unique = residualize_primary_relation(
        source,
        target,
        primary,
        controls,
    )
    support = dyadic_signal_support(source, target, residual)

    simulation = contract["detectability_surface"]
    fixed_kind = "environment" if args.kind in {"environment", "current"} else "historical"
    fixed_map = simulation["fixed_control_effects"][fixed_kind]
    fixed = np.zeros(len(predictor_names), dtype=float)
    index = {name: i for i, name in enumerate(predictor_names)}
    for name, value in fixed_map.items():
        fixed[index[name]] = float(value)

    worlds = (
        int(args.worlds_per_cell)
        if args.worlds_per_cell is not None
        else int(simulation["worlds_per_cell"])
    )
    surface = detectability_surface(
        source,
        target,
        predictors,
        effects=simulation["standardized_effect_grid"],
        private_amplitudes=simulation["private_amplitudes"],
        primary_index=0,
        alpha=float(simulation["intended_inference_alpha"]),
        worlds_per_cell=worlds,
        source_intercept_sd=float(simulation["source_intercept_sd"]),
        target_intercept_sd=float(simulation["target_intercept_sd"]),
        dyad_noise_sd=float(simulation["dyad_noise_sd"]),
        fixed_coefficients=fixed,
        source_random_slope_index=index["z_geographic_coverage"],
        target_random_slope_index=index["centered_same_order"],
        random_slope_sd_per_amplitude=float(
            simulation["random_slope_sd_per_amplitude"]
        ),
        response_transform=str(simulation["response_transform"]),
        seed_namespace=(
            f"{simulation['seed_namespace']}|{args.kind}|{args.panel}|"
            f"{sha256_path(args.design)}"
        ),
        block_size=int(simulation["block_size"]),
    )

    repeatability = None
    if args.relation_replicates is not None:
        relation_matrix = np.load(
            args.relation_replicates,
            allow_pickle=False,
        )
        repeatability = asdict(relation_repeatability(relation_matrix))

    payload = {
        "schema": "ttf_q_characterization_result_v0.1",
        "status": "CHARACTERIZED_RESPONSE_BLIND_PREDICTOR_DOMAIN",
        "kind": args.kind,
        "panel": args.panel,
        "inputs": {
            "design_sha256": sha256_path(args.design),
            "contract_sha256": sha256_path(args.contract),
            "relation_replicates_sha256": (
                sha256_path(args.relation_replicates)
                if args.relation_replicates is not None
                else None
            ),
        },
        "predictors_in_order": predictor_names,
        "relation_repeatability": repeatability,
        "nonredundancy": asdict(unique),
        "dyadic_signal_support": asdict(support),
        "detectability_surface": surface,
        "minimum_detectable_effect_by_private_amplitude": (
            minimum_detectable_effects(
                surface,
                target_power=float(simulation["target_power"]),
            )
        ),
        "calibrated_evaluable_envelope": calibrated_detectability_envelope(
            surface,
            target_power=float(simulation["target_power"]),
            type1_wilson_upper_max=0.05,
        ),
        "interpretation_policy": {
            "binary_pass_fail_forbidden": True,
            "biological_effect_tested": False,
            "genetic_response_opened": False,
            "meaning": (
                "Describes the response-blind domain in which the focal pairwise "
                "ecological relation is measured distinctly enough and with enough "
                "design support to sustain intended inference. It is not a result "
                "about whether the biological transferability effect exists."
            ),
        },
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n"
    )
    print(
        json.dumps(
            {
                "status": payload["status"],
                "kind": args.kind,
                "panel": args.panel,
                "nonredundancy": payload["nonredundancy"],
                "dyadic_signal_support": payload["dyadic_signal_support"],
                "minimum_detectable_effect_by_private_amplitude": payload[
                    "minimum_detectable_effect_by_private_amplitude"
                ],
            },
            sort_keys=True,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
