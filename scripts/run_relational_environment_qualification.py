#!/usr/bin/env python3
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np

from ttf.relational_dyadic import (
    batch_primary_test,
    prepare_dyadic_regression,
    wilson_interval,
)


def sha256_path(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def zscore(x):
    x = np.asarray(x, float)
    sd = float(x.std())
    if not np.isfinite(sd) or sd <= np.finfo(float).eps:
        raise ValueError("constant predictor")
    return (x - x.mean()) / sd


def remap(x):
    vals = sorted(map(int, np.unique(x)))
    lookup = {v: i for i, v in enumerate(vals)}
    return np.asarray([lookup[int(v)] for v in x], np.int64)


def seed_for(namespace, design_sha, cell):
    return int(
        hashlib.sha256(
            f"{namespace}|{design_sha}|{cell}".encode()
        ).hexdigest()[:16],
        16,
    )


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--opportunity-design", type=Path, required=True)
    ap.add_argument("--opportunity-summary", type=Path, required=True)
    ap.add_argument("--rule", type=Path, required=True)
    ap.add_argument("--output", type=Path, required=True)
    ap.add_argument("--block-size", type=int, default=100)
    args = ap.parse_args()

    rule = json.loads(args.rule.read_text())
    summary = json.loads(args.opportunity_summary.read_text())
    if rule.get("schema") != "ttf_relational_environment_qualification_rule_v0.2":
        raise RuntimeError("bad rule")
    if summary.get("status") != "PASS_TO_DEVELOPMENT_SYNTHETIC_QUALIFICATION":
        raise RuntimeError("development opportunity gate not passed")
    if any(bool(v) for v in summary["response_firewall"].values()):
        raise RuntimeError("Study B response firewall open")

    data = np.load(args.opportunity_design, allow_pickle=False)
    s = np.asarray(data["development_source_index"], np.int64)
    t = np.asarray(data["development_target_index"], np.int64)
    r = np.asarray(data["development_R_env"], float)
    g = np.asarray(data["development_coverage"], float)
    sc = np.asarray(data["development_same_class"], float)
    so = np.asarray(data["development_same_order"], float)
    sf = np.asarray(data["development_same_family"], float)
    lr = np.asarray(data["development_locality_count_ratio"], float)

    x = np.column_stack((
        zscore(r),
        zscore(g),
        sc - sc.mean(),
        so - so.mean(),
        sf - sf.mean(),
        zscore(np.abs(np.log(lr))),
    ))
    predictor_names = list(map(str, rule["inference"]["predictors_in_order"]))
    if len(predictor_names) != x.shape[1]:
        raise RuntimeError("predictor contract width drift")
    predictor_index = {name: i for i, name in enumerate(predictor_names)}
    primary_index = int(rule["inference"]["primary_index"])
    if primary_index != predictor_index.get("z_R_env"):
        raise RuntimeError("primary predictor contract drift")

    source_group = remap(s)
    target_group = remap(t)
    prepared = prepare_dyadic_regression(
        source_group, target_group, x, primary_index=primary_index
    )

    synth = rule["synthetic_worlds"]
    gate = rule["gate_numeric"]
    condition_max = float(gate["predictor_condition_number_max_exclusive"])
    if prepared.condition_number >= condition_max:
        raise RuntimeError("condition-number gate failed")

    worlds = int(synth["worlds_per_cell"])
    design_sha = sha256_path(args.opportunity_design)
    alpha = float(gate["p_value_cutoff"])
    if not np.isclose(alpha, float(rule["inference"]["alpha"])):
        raise RuntimeError("qualification alpha drift")
    type1_upper = float(gate["private_type1_wilson95_upper_max"])
    power_lower = float(gate["relational_power_wilson95_lower_min"])

    source_intercept_sd = float(synth["source_intercept_sd"])
    target_intercept_sd = float(synth["target_intercept_sd"])
    noise_sd = float(synth["dyad_noise_sd"])

    fixed = np.zeros(x.shape[1], dtype=float)
    for name, value in synth["fixed_control_effects"].items():
        if name not in predictor_index:
            raise RuntimeError(f"unknown frozen fixed-control predictor: {name}")
        if predictor_index[name] == primary_index:
            raise RuntimeError("primary predictor cannot be a fixed nuisance control")
        fixed[predictor_index[name]] = float(value)

    private = synth["private_structure_numeric"]
    private_sd_per_A = float(private["slope_sd_per_A"])
    source_private_index = predictor_index[str(private["source_predictor"])]
    target_private_index = predictor_index[str(private["target_predictor"])]
    if primary_index in (source_private_index, target_private_index):
        raise RuntimeError("private nuisance slope cannot target primary predictor")

    ns = prepared.absorber.n_source
    nt = prepared.absorber.n_target
    n = len(s)
    block = int(args.block_size)
    cells = list(synth["private_null_cells"]) + [
        dict(synth["relational_positive_cell"])
    ]

    out = {}
    fixed_term = x @ fixed
    for cell in cells:
        name = str(cell["name"])
        A = float(cell["A"])
        beta = float(cell["beta_R_latent"])
        rng = np.random.default_rng(
            seed_for(synth["seed_namespace"], design_sha, name)
        )
        rejected = 0
        coefs = []
        ses = []
        completed = 0
        while completed < worlds:
            w = min(block, worlds - completed)
            source_intercept = rng.normal(
                0, source_intercept_sd, size=(ns, w)
            )
            target_intercept = rng.normal(
                0, target_intercept_sd, size=(nt, w)
            )
            source_private = rng.normal(
                0, private_sd_per_A * A, size=(ns, w)
            )
            target_private = rng.normal(
                0, private_sd_per_A * A, size=(nt, w)
            )
            noise = rng.normal(0, noise_sd, size=(n, w))
            latent = (
                source_intercept[source_group]
                + target_intercept[target_group]
                + fixed_term[:, None]
                + source_private[source_group]
                * x[:, source_private_index, None]
                + target_private[target_group]
                * x[:, target_private_index, None]
                + noise
                + beta * x[:, primary_index, None]
            )
            fit = batch_primary_test(prepared, np.tanh(latent))
            if not (
                np.isfinite(fit.coefficient).all()
                and np.isfinite(fit.standard_error).all()
                and np.isfinite(fit.p_value_one_sided).all()
            ):
                raise RuntimeError(f"non-finite world in {name}")
            rejected += int(
                np.count_nonzero(fit.p_value_one_sided <= alpha)
            )
            coefs.extend(map(float, fit.coefficient))
            ses.extend(map(float, fit.standard_error))
            completed += w

        lo, hi = wilson_interval(rejected, worlds)
        out[name] = {
            "A": A,
            "beta_R_latent": beta,
            "worlds": worlds,
            "alpha": alpha,
            "rejections_p_le_alpha": rejected,
            "rejection_rate": rejected / worlds,
            "wilson95_lower": lo,
            "wilson95_upper": hi,
            "coefficient_mean": float(np.mean(coefs)),
            "coefficient_sd": float(np.std(coefs)),
            "standard_error_median": float(np.median(ses)),
            "master_seed_uint64": seed_for(
                synth["seed_namespace"], design_sha, name
            ),
        }

    type1 = all(
        out[c["name"]]["wilson95_upper"] <= type1_upper
        for c in synth["private_null_cells"]
    )
    pos = synth["relational_positive_cell"]["name"]
    power = out[pos]["wilson95_lower"] >= power_lower
    status = (
        "PASS_TO_CONFIRMATORY_CHARACTER_MASK_PREPARATION"
        if type1 and power
        else "NOT_EVALUABLE_ENVIRONMENT_SYNTHETIC_QUALIFICATION"
    )
    payload = {
        "schema": "ttf_relational_environment_qualification_result_v0.2",
        "status": status,
        "rule_sha256": sha256_path(args.rule),
        "opportunity_design_sha256": design_sha,
        "alpha": alpha,
        "geometry": {
            "dyads": n,
            "source_clusters": ns,
            "target_clusters": nt,
            "predictor_condition_number": prepared.condition_number,
        },
        "world_contract": {
            "source_intercept_sd": source_intercept_sd,
            "target_intercept_sd": target_intercept_sd,
            "dyad_noise_sd": noise_sd,
            "fixed_control_effects": synth["fixed_control_effects"],
            "private_structure_numeric": private,
        },
        "cells": out,
        "gates": {
            "predictor_condition_number_max_exclusive": condition_max,
            "private_type1_wilson95_upper_max": type1_upper,
            "relational_power_wilson95_lower_min": power_lower,
            "private_type1_pass": type1,
            "relational_power_pass": power,
            "overall_pass": bool(type1 and power),
        },
        "response_firewall": {
            "Study_B_sequence_identity_opened": False,
            "Study_B_pairwise_genetic_distances_opened": False,
            "Study_B_T_st_computed": False,
            "Study_B_beta_R_computed": False,
        },
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n"
    )
    print(json.dumps({
        "status": status,
        "gates": payload["gates"],
        "geometry": payload["geometry"],
        "cells": out,
    }, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
