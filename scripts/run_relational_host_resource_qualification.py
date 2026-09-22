#!/usr/bin/env python3
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import numpy as np

from ttf.relational_dyadic import batch_primary_test, prepare_dyadic_regression, wilson_interval


def sha256_path(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def zscore(x: np.ndarray) -> np.ndarray:
    x = np.asarray(x, dtype=float)
    sd = float(x.std())
    if not np.isfinite(sd) or sd <= np.finfo(float).eps:
        raise ValueError("cannot z-score constant/non-finite predictor")
    return (x - float(x.mean())) / sd


def seed_for(namespace: str, design_sha: str, cell: str) -> int:
    digest = hashlib.sha256(f"{namespace}|{design_sha}|{cell}".encode()).hexdigest()
    return int(digest[:16], 16)


def remap(values: np.ndarray) -> tuple[np.ndarray, list[int]]:
    unique = sorted(map(int, np.unique(values)))
    mapping = {value: index for index, value in enumerate(unique)}
    return np.asarray([mapping[int(v)] for v in values], dtype=np.int64), unique


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--design-npz", type=Path, required=True)
    ap.add_argument("--audit", type=Path, required=True)
    ap.add_argument("--rule", type=Path, required=True)
    ap.add_argument("--output", type=Path, required=True)
    ap.add_argument("--block-size", type=int, default=100)
    args = ap.parse_args()

    rule = json.loads(args.rule.read_text())
    if rule.get("schema") != "ttf_relational_host_resource_qualification_rule_v0.1":
        raise RuntimeError("unexpected qualification rule")
    if sha256_path(args.design_npz) != rule["design_source"]["sha256"]:
        raise RuntimeError("host-resource design SHA256 drift")
    audit = json.loads(args.audit.read_text())
    if audit.get("schema") != "ttf_relational_host_resource_response_blind_audit_v0.1":
        raise RuntimeError("unexpected response-blind audit")
    if any(bool(v) for v in audit["response_firewall"].values()):
        raise RuntimeError("audit response firewall is open")

    data = np.load(args.design_npz, allow_pickle=False)
    species = np.asarray(data["species_order"]).astype(str)
    family = np.asarray(data["family"]).astype(str)
    source0 = np.asarray(data["source_index"], dtype=np.int64)
    target0 = np.asarray(data["target_index"], dtype=np.int64)
    dropped = {str(row["species"]) for row in audit["freshness"]["dropped_species"]}
    drop_flag = np.asarray([name in dropped for name in species], dtype=bool)
    mask = ~(drop_flag[source0] | drop_flag[target0])
    source0 = source0[mask]
    target0 = target0[mask]

    if int(mask.sum()) != int(rule["design_source"]["supported_directed_pairs"]):
        raise RuntimeError("supported pair count drift after frozen freshness drop")

    source, source_levels = remap(source0)
    target, target_levels = remap(target0)
    if len(source_levels) != int(rule["design_source"]["supported_source_clusters"]):
        raise RuntimeError("source cluster count drift")
    if len(target_levels) != int(rule["design_source"]["supported_target_clusters"]):
        raise RuntimeError("target cluster count drift")

    r = zscore(np.asarray(data["host_resource_jaccard"], dtype=float)[mask])
    g = zscore(np.asarray(data["coverage"], dtype=float)[mask])
    d = zscore(np.log1p(np.asarray(data["centroid_distance_km"], dtype=float)[mask]))
    m = zscore(np.abs(np.log(np.asarray(data["locality_count_ratio"], dtype=float)[mask])))
    f = (family[source0] == family[target0]).astype(float)
    f = f - float(f.mean())
    x = np.column_stack((r, g, d, m, f))
    prepared = prepare_dyadic_regression(source, target, x, primary_index=0)
    if prepared.condition_number >= 1e4:
        raise RuntimeError("response-blind predictor condition number gate failed")

    worlds_per_cell = int(rule["synthetic_worlds"]["worlds_per_cell"])
    block_size = int(args.block_size)
    if worlds_per_cell < 1 or block_size < 1:
        raise ValueError("world and block sizes must be positive")
    namespace = str(rule["synthetic_worlds"]["seed_namespace"])
    design_sha = str(rule["design_source"]["sha256"])
    ns = prepared.absorber.n_source
    nt = prepared.absorber.n_target
    n = len(source)

    cells = list(rule["synthetic_worlds"]["private_null_cells"]) + [
        dict(rule["synthetic_worlds"]["relational_positive_cell"])
    ]
    results = {}
    for cell in cells:
        name = str(cell["name"])
        amplitude = float(cell["A"])
        beta_r = float(cell["beta_R_latent"])
        rng = np.random.default_rng(seed_for(namespace, design_sha, name))
        rejected = 0
        coefficients = []
        standard_errors = []
        completed = 0
        while completed < worlds_per_cell:
            width = min(block_size, worlds_per_cell - completed)
            source_intercept = rng.normal(0.0, 0.18, size=(ns, width))
            target_intercept = rng.normal(0.0, 0.18, size=(nt, width))
            source_private = rng.normal(0.0, 0.12 * amplitude, size=(ns, width))
            target_private = rng.normal(0.0, 0.12 * amplitude, size=(nt, width))
            dyad_noise = rng.normal(0.0, 0.32, size=(n, width))
            latent = (
                source_intercept[source]
                + target_intercept[target]
                + 0.10 * x[:, 1, None]
                - 0.08 * x[:, 2, None]
                + 0.05 * x[:, 3, None]
                + 0.04 * x[:, 4, None]
                + source_private[source] * x[:, 1, None]
                + target_private[target] * x[:, 2, None]
                + dyad_noise
                + beta_r * x[:, 0, None]
            )
            response = np.tanh(latent)
            fitted = batch_primary_test(prepared, response)
            if not (
                np.isfinite(fitted.coefficient).all()
                and np.isfinite(fitted.standard_error).all()
                and np.isfinite(fitted.p_value_one_sided).all()
            ):
                raise RuntimeError(f"non-finite formal result in {name}")
            rejected += int(np.count_nonzero(fitted.p_value_one_sided <= 0.05))
            coefficients.extend(map(float, fitted.coefficient))
            standard_errors.extend(map(float, fitted.standard_error))
            completed += width

        lower, upper = wilson_interval(rejected, worlds_per_cell)
        results[name] = {
            "A": amplitude,
            "beta_R_latent": beta_r,
            "worlds": worlds_per_cell,
            "rejections_p_le_0p05": rejected,
            "rejection_rate": rejected / worlds_per_cell,
            "wilson95_lower": lower,
            "wilson95_upper": upper,
            "coefficient_mean": float(np.mean(coefficients)),
            "coefficient_sd": float(np.std(coefficients)),
            "standard_error_median": float(np.median(standard_errors)),
            "master_seed_uint64": seed_for(namespace, design_sha, name),
        }

    null_pass = all(
        results[str(cell["name"])]["wilson95_upper"] <= 0.10
        for cell in rule["synthetic_worlds"]["private_null_cells"]
    )
    positive_name = str(rule["synthetic_worlds"]["relational_positive_cell"]["name"])
    power_pass = results[positive_name]["wilson95_lower"] >= 0.80
    overall = "PASS_TO_CHARACTER_MASK_GATE" if null_pass and power_pass else "NOT_EVALUABLE_SYNTHETIC_QUALIFICATION"

    payload = {
        "schema": "ttf_relational_host_resource_qualification_result_v0.1",
        "status": overall,
        "rule_sha256": sha256_path(args.rule),
        "audit_sha256": sha256_path(args.audit),
        "design_npz_sha256": sha256_path(args.design_npz),
        "geometry": {
            "dyads": n,
            "source_clusters": ns,
            "target_clusters": nt,
            "predictor_condition_number": prepared.condition_number,
        },
        "cells": results,
        "gates": {
            "private_type1_pass": null_pass,
            "relational_power_pass": power_pass,
            "overall_pass": bool(null_pass and power_pass),
        },
        "next_step": (
            "Open only the canonical-valid/noncanonical character mask for the frozen 487-species host panel, then rerun this exact qualification on the survivor geometry before any nucleotide identity or T_st."
            if null_pass and power_pass
            else "STOP. Do not open character masks, nucleotide identity, pairwise genetic distances, or T_st."
        ),
        "response_firewall": {
            "sequence_identity_opened": False,
            "pairwise_genetic_distances_opened": False,
            "T_st_computed": False,
            "beta_R_computed": False,
        },
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"status": overall, "gates": payload["gates"], "cells": results}, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
