#!/usr/bin/env python3
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np

from ttf.precision import wilson_interval


RULE_SCHEMA = "ttf_genetic_conditional_order_contrast_qualification_rule_v0.2"
SHARD_SCHEMA = "ttf_genetic_conditional_order_contrast_qualification_shard_v0.2"


def sha256_path(path: Path) -> str:
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def load_json(path: Path, schema: str) -> dict:
    payload = json.loads(Path(path).read_text())
    if payload.get("schema") != schema:
        raise RuntimeError(f"unexpected schema for {path}: {payload.get('schema')!r}")
    return payload


def expected_cells(rule: dict) -> dict[str, dict]:
    out = {
        str(row["cell"]): dict(row)
        for row in rule["synthetic_worlds"]["private_null_cells"]
    }
    positive = dict(rule["synthetic_worlds"]["positive_control"])
    out[str(positive["cell"])] = positive
    return out


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--rule", type=Path, required=True)
    parser.add_argument("--shard", type=Path, action="append", required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    rule = load_json(args.rule, RULE_SCHEMA)
    actual_rule_sha = sha256_path(args.rule)
    cells = expected_cells(rule)
    by_cell: dict[str, dict[int, dict]] = {name: {} for name in cells}
    metadata: dict[str, tuple[str, str, int]] = {}

    for path in args.shard:
        shard = load_json(path, SHARD_SCHEMA)
        if shard.get("status") != "SYNTHETIC_ORDER_CONTRAST_QUALIFICATION_SHARD_COMPLETE":
            raise RuntimeError(f"incomplete shard: {path}")
        if shard["cell"] not in cells:
            raise RuntimeError(f"unexpected cell in {path}: {shard['cell']}")
        if any(value is not False for value in shard["outcome_firewall"].values()):
            raise RuntimeError(f"open empirical outcome firewall in {path}")
        if str(shard["rule_sha256"]) != actual_rule_sha:
            raise RuntimeError("qualification shard rule hash does not match supplied frozen rule")
        meta = (
            str(shard["rule_sha256"]),
            str(shard["geometry_csv_sha256"]),
            int(shard["n_eval_species"]),
        )
        previous = metadata.setdefault(
            str(shard["geometry_manifest_sha256"]),
            meta,
        )
        if previous != meta:
            raise RuntimeError("qualification shard metadata drift")
        for row in shard["rows"]:
            index = int(row["absolute_replicate_index"])
            if index in by_cell[shard["cell"]]:
                raise RuntimeError(
                    f"duplicate replicate {shard['cell']}:{index}"
                )
            if not np.isfinite(float(row["statistic"])):
                raise RuntimeError("non-finite qualification statistic")
            p = float(row["p_value"])
            if not 0.0 < p <= 1.0:
                raise RuntimeError("invalid qualification p-value")
            if bool(row["reject"]) is not bool(
                p <= float(rule["frozen_estimator"]["alpha"])
            ):
                raise RuntimeError("shard rejection label drift")
            by_cell[shard["cell"]][index] = dict(row)

    if len(metadata) != 1:
        raise RuntimeError("qualification shards do not share one geometry manifest")
    geometry_manifest_sha, (
        rule_sha,
        geometry_csv_sha,
        n_eval_species,
    ) = next(iter(metadata.items()))
    if rule_sha != actual_rule_sha:
        raise RuntimeError("qualification shard rule hash does not match supplied frozen rule")

    results = {}
    for name, spec in cells.items():
        expected_n = int(spec["worlds"])
        observed = by_cell[name]
        expected_indices = set(range(expected_n))
        if set(observed) != expected_indices:
            missing = sorted(expected_indices - set(observed))
            extra = sorted(set(observed) - expected_indices)
            raise RuntimeError(
                f"incomplete cell {name}: missing={missing[:10]}, extra={extra[:10]}"
            )
        rows = [observed[i] for i in range(expected_n)]
        rejections = int(sum(bool(row["reject"]) for row in rows))
        interval = wilson_interval(rejections, expected_n)
        results[name] = {
            "worlds": expected_n,
            "rejections": rejections,
            "rejection_rate": float(rejections / expected_n),
            "wilson95_low": interval.low,
            "wilson95_high": interval.high,
            "mean_statistic": float(np.mean([float(row["statistic"]) for row in rows])),
            "median_p_value": float(np.median([float(row["p_value"]) for row in rows])),
        }

    private_names = [
        str(row["cell"]) for row in rule["synthetic_worlds"]["private_null_cells"]
    ]
    positive_name = str(rule["synthetic_worlds"]["positive_control"]["cell"])
    type1_ceiling = float(rule["qualification"]["type1_wilson95_upper_ceiling"])
    power_floor = float(rule["qualification"]["power_wilson95_lower_floor"])
    type1_pass = all(
        results[name]["wilson95_high"] <= type1_ceiling
        for name in private_names
    )
    power_pass = results[positive_name]["wilson95_low"] >= power_floor
    passed = bool(type1_pass and power_pass)

    payload = {
        "schema": "ttf_genetic_conditional_order_contrast_qualification_v0.2",
        "status": "PASS" if passed else "NOT_EVALUABLE_CONDITIONAL_ORDER_CONTRAST_V02",
        "passed": passed,
        "rule_sha256": rule_sha,
        "geometry_manifest_sha256": geometry_manifest_sha,
        "geometry_csv_sha256": geometry_csv_sha,
        "n_eval_species": n_eval_species,
        "type1_gate": {
            "pass": bool(type1_pass),
            "wilson95_upper_ceiling": type1_ceiling,
            "max_observed_wilson95_upper": max(
                results[name]["wilson95_high"] for name in private_names
            ),
        },
        "power_gate": {
            "pass": bool(power_pass),
            "wilson95_lower_floor": power_floor,
            "observed_wilson95_lower": results[positive_name]["wilson95_low"],
        },
        "cells": results,
        "confirmatory_identity_opening_eligible": False,
        "next_step": (
            "v0.2 contrast qualification PASS only; confirmatory character-mask admissibility and exact survivor-geometry requalification are still required before any nucleotide identity opening"
            if passed
            else "STOP; do not open confirmatory nucleotide identity"
        ),
        "outcome_firewall": {
            "development_sequence_identity_opened": False,
            "confirmatory_sequence_identity_opened": False,
            "confirmatory_pairwise_genetic_distances_opened": False,
            "conditional_empirical_result_opened": False,
        },
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    print(
        json.dumps(
            {
                "status": payload["status"],
                "type1_pass": type1_pass,
                "power_pass": power_pass,
                "max_type1_upper": payload["type1_gate"]["max_observed_wilson95_upper"],
                "power_lower": payload["power_gate"]["observed_wilson95_lower"],
            },
            sort_keys=True,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
