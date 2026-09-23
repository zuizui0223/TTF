#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from collections import Counter
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


def species_digest(names: list[str]) -> str:
    return hashlib.sha256(("\n".join(sorted(names)) + "\n").encode()).hexdigest()


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--design-npz", type=Path, required=True)
    ap.add_argument("--survivors", type=Path, required=True)
    ap.add_argument("--repair", type=Path, required=True)
    ap.add_argument("--rule", type=Path, required=True)
    ap.add_argument("--output", type=Path, required=True)
    ap.add_argument("--block-size", type=int, default=100)
    args = ap.parse_args()

    repair = json.loads(args.repair.read_text())
    if repair.get("schema") != "ttf_relational_host_resource_mask_parser_repair_v0.1":
        raise RuntimeError("unexpected repair contract")
    if repair["response_state_at_repair_freeze"] != {
        "T_st_opened": False,
        "beta_R_opened": False,
        "relational_decision_opened": False,
    }:
        raise RuntimeError("repair contract was frozen after relational response opening")

    rule = json.loads(args.rule.read_text())
    if rule.get("schema") != "ttf_relational_host_resource_qualification_rule_v0.1":
        raise RuntimeError("unexpected qualification rule")
    if sha256_path(args.design_npz) != rule["design_source"]["sha256"]:
        raise RuntimeError("host-resource design SHA256 drift")

    rows = list(csv.DictReader(args.survivors.open(encoding="utf-8")))
    original_survivors = {str(row["species"]).strip() for row in rows}
    excluded = {
        str(row["species"])
        for row in repair["response_blind_header_audit"]["affected_species"]
    }
    survivors = original_survivors - excluded
    expected_geometry = repair["repaired_survivor_geometry_contract"]
    if len(survivors) != int(expected_geometry["repaired_mask_survivors"]):
        raise RuntimeError("repaired survivor count drift")
    if species_digest(list(survivors)) != expected_geometry["repaired_survivor_species_digest_sha256"]:
        raise RuntimeError("repaired survivor digest drift")

    data = np.load(args.design_npz, allow_pickle=False)
    species = np.asarray(data["species_order"]).astype(str)
    family = np.asarray(data["family"]).astype(str)
    source0_all = np.asarray(data["source_index"], dtype=np.int64)
    target0_all = np.asarray(data["target_index"], dtype=np.int64)

    pair_mask = np.asarray(
        [species[s] in survivors and species[t] in survivors for s, t in zip(source0_all, target0_all)],
        dtype=bool,
    )
    target_counts = Counter(map(int, target0_all[pair_mask]))
    eligible_targets = {
        target for target, count in target_counts.items()
        if count >= int(repair["repaired_survivor_geometry_contract"]["minimum_sources_per_target"])
    }
    pair_mask &= np.asarray([int(t) in eligible_targets for t in target0_all], dtype=bool)

    source0 = source0_all[pair_mask]
    target0 = target0_all[pair_mask]
    source_names = sorted(set(map(str, species[source0])))
    target_names = sorted(set(map(str, species[target0])))
    union_names = sorted(set(source_names) | set(target_names))

    if int(pair_mask.sum()) != int(expected_geometry["repaired_supported_directed_pairs"]):
        raise RuntimeError("repaired directed-pair count drift")
    if len(source_names) != int(expected_geometry["repaired_supported_source_clusters"]):
        raise RuntimeError("repaired source cluster count drift")
    if len(target_names) != int(expected_geometry["repaired_supported_target_clusters"]):
        raise RuntimeError("repaired target cluster count drift")
    if species_digest(source_names) != expected_geometry["source_digest_sha256"]:
        raise RuntimeError("repaired source digest drift")
    if species_digest(target_names) != expected_geometry["target_digest_sha256"]:
        raise RuntimeError("repaired target digest drift")
    if species_digest(union_names) != expected_geometry["empirical_union_digest_sha256"]:
        raise RuntimeError("repaired empirical union digest drift")

    source, source_levels = remap(source0)
    target, target_levels = remap(target0)

    r = zscore(np.asarray(data["host_resource_jaccard"], dtype=float)[pair_mask])
    g = zscore(np.asarray(data["coverage"], dtype=float)[pair_mask])
    d = zscore(np.log1p(np.asarray(data["centroid_distance_km"], dtype=float)[pair_mask]))
    m = zscore(np.abs(np.log(np.asarray(data["locality_count_ratio"], dtype=float)[pair_mask])))
    f = (family[source0] == family[target0]).astype(float)
    f = f - float(f.mean())
    x = np.column_stack((r, g, d, m, f))

    prepared = prepare_dyadic_regression(source, target, x, primary_index=0)
    if prepared.condition_number >= 1e4:
        raise RuntimeError("response-blind predictor condition-number gate failed")

    worlds_per_cell = int(rule["synthetic_worlds"]["worlds_per_cell"])
    namespace = str(rule["synthetic_worlds"]["seed_namespace"])
    design_sha = str(rule["design_source"]["sha256"])
    block_size = int(args.block_size)
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
            fitted = batch_primary_test(prepared, np.tanh(latent))
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
    overall = (
        "PASS_TO_REPAIRED_EMPIRICAL_AUTHORIZATION"
        if null_pass and power_pass
        else "NOT_EVALUABLE_REPAIRED_SURVIVOR_GEOMETRY"
    )

    payload = {
        "schema": "ttf_relational_host_resource_repaired_survivor_requalification_v0.1",
        "status": overall,
        "repair_sha256": sha256_path(args.repair),
        "rule_sha256": sha256_path(args.rule),
        "design_npz_sha256": sha256_path(args.design_npz),
        "geometry": {
            "mask_survivors": len(survivors),
            "dyads": n,
            "source_clusters": ns,
            "target_clusters": nt,
            "empirical_union_species": len(union_names),
            "predictor_condition_number": prepared.condition_number,
            "source_digest_sha256": species_digest(source_names),
            "target_digest_sha256": species_digest(target_names),
            "empirical_union_digest_sha256": species_digest(union_names),
        },
        "cells": results,
        "gates": {
            "private_type1_pass": null_pass,
            "relational_power_pass": power_pass,
            "overall_pass": bool(null_pass and power_pass),
        },
        "response_state": {
            "T_st_computed": False,
            "beta_R_computed": False,
            "relational_decision_opened": False,
        },
        "next_step": (
            "Issue a new one-continuation empirical authorization bound to the repaired 420-species mask-survivor contract."
            if null_pass and power_pass
            else "STOP. Do not continue empirical relational scoring."
        ),
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    print(json.dumps({
        "status": overall,
        "geometry": payload["geometry"],
        "gates": payload["gates"],
    }, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
