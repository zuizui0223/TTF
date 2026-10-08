#!/usr/bin/env python3
"""Response-blind, climate-valid HOST cloud support census, before genetic opening.

This only checks whether *frozen* larval host-resource climate clouds can be
constructed from exact WCVP/HOSTS/WGSRPD3/CHELSA 0-BP sources. Passing this
partial external-data gate never authorizes an empirical genetic result.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from collections import defaultdict
from pathlib import Path

import numpy as np

try:
    from scripts.build_historical_host_memory_predictor import (
        FROZEN_INPUT_SHA256, norm_id, sample_matrix, sha256_path,
        verify_host_pair_identity,
    )
except ModuleNotFoundError:
    from build_historical_host_memory_predictor import (
        FROZEN_INPUT_SHA256, norm_id, sample_matrix, sha256_path,
        verify_host_pair_identity,
    )


VARIABLE_ORDER = ("bio01", "bio07", "bio12", "bio15")


def census_finite_host_support(candidates, native_rows, support_rows, finite):
    """Count climate-finite WGSRPD3 sample points in each frozen host union."""
    keep = np.asarray(finite, dtype=bool)
    if keep.ndim != 1 or len(keep) != len(support_rows):
        raise RuntimeError("support point/finite climate mask length mismatch")
    support_counts = defaultdict(int)
    for row, is_finite in zip(support_rows, keep):
        area = str(row["area_code_l3"]).strip()
        if not area:
            raise RuntimeError("empty WGSRPD3 support region")
        support_counts[area] += int(is_finite)

    units_by_host = defaultdict(set)
    for row in native_rows:
        area = str(row["area_code_l3"]).strip()
        if not area:
            raise RuntimeError("empty native WGSRPD3 area")
        units_by_host[norm_id(row["accepted_plant_name_id"])].add(area)

    counts = []
    for row in candidates:
        species = row["species"].strip()
        host_ids = [norm_id(x) for x in row["accepted_host_ids"].split(";")]
        units = set().union(*(units_by_host.get(h, set()) for h in host_ids))
        n = sum(support_counts.get(area, 0) for area in units)
        counts.append({"species": species, "host_climate_points": int(n),
                       "has_minimum_ten": n >= 10})
    return counts


def write_receipt(candidates_path, pairs_path, native_path, support_path,
                  raster_paths, counts, support_finite_count,
                  rule_path, output_json, output_csv):
    rule = json.loads(rule_path.read_text(encoding="utf-8"))
    if rule.get("schema") != "ttf_historical_host_memory_predictor_materialization_rule_v0.1":
        raise RuntimeError("wrong frozen host-memory materialization rule")
    if any(rule.get("response_firewall", {}).values()):
        raise RuntimeError("frozen response firewall not fully closed")
    minimum = int(rule["host_cloud"]["minimum_Q_host_points"])
    if minimum != 10:
        raise RuntimeError("frozen minimum host climate cloud changed")
    if len(counts) != 642:
        raise RuntimeError("expected 642 fixed species, no backfill")
    eligible = sum(x["host_climate_points"] >= minimum for x in counts)
    dest = Path(output_csv)
    dest.parent.mkdir(parents=True, exist_ok=True)
    with dest.open("w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=["species", "host_climate_points",
                                               "has_minimum_ten"], lineterminator="\n")
        writer.writeheader()
        writer.writerows(counts)

    receipt = {
        "schema": "ttf_historical_host_memory_finite_host_cloud_census_v0.1",
        "status": ("PASS_PARTIAL_HOST_CLOUD_CLIMATE_SUPPORT" if eligible >= 500
                   else "NOT_EVALUABLE_HOST_CLOUD_CLIMATE_SUPPORT"),
        "candidate_species": 642,
        "host_cloud_valid_species": eligible,
        "host_cloud_invalid_species": 642 - eligible,
        "minimum_host_points": minimum,
        "minimum_species_for_full_predictor": 500,
        "finite_wgsrpd3_support_points": int(support_finite_count),
        "host_cloud_point_min": int(min(x["host_climate_points"] for x in counts)),
        "host_cloud_point_max": int(max(x["host_climate_points"] for x in counts)),
        "input_sha256": {
            "candidates": sha256_path(candidates_path),
            "host_pairs": sha256_path(pairs_path),
            "native_units": sha256_path(native_path),
            "wgsrpd_support": sha256_path(support_path),
            "materialization_rule": sha256_path(rule_path),
            **{p.name: sha256_path(p) for p in raster_paths},
        },
        "species_support_csv_sha256": sha256_path(dest),
        "interpretation": "Partial external-data qualification only. No self-history, edge predictors, synthetic qualifications, or genetic responses computed.",
        "full_predictor_qualification_passed": False,
        "genetic_response_opened": False,
    }
    output_json = Path(output_json)
    output_json.parent.mkdir(parents=True, exist_ok=True)
    output_json.write_text(json.dumps(receipt, indent=2, sort_keys=True) + "\n",
                           encoding="utf-8")
    print(json.dumps({"status": receipt["status"],
                      "valid_species": eligible,
                      "finite_support": int(support_finite_count)}, sort_keys=True))
    return receipt


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--candidates", type=Path, required=True)
    ap.add_argument("--host-pairs", type=Path, required=True)
    ap.add_argument("--native-units", type=Path, required=True)
    ap.add_argument("--wgsrpd-support", type=Path, required=True)
    ap.add_argument("--rule", type=Path, required=True)
    ap.add_argument("--current-asset", type=Path, required=True, action="append")
    ap.add_argument("--output-receipt", type=Path, required=True)
    ap.add_argument("--output-species", type=Path, required=True)
    args = ap.parse_args()

    for key, path in (("candidates", args.candidates),
                      ("host_pairs", args.host_pairs),
                      ("native_units", args.native_units),
                      ("wgsrpd_support", args.wgsrpd_support)):
        if sha256_path(path) != FROZEN_INPUT_SHA256[key]:
            raise RuntimeError(f"exact frozen {key} source mismatch")

    with args.candidates.open(newline="", encoding="utf-8") as f:
        candidates = list(csv.DictReader(f))
    with args.host_pairs.open(newline="", encoding="utf-8") as f:
        pair_rows = list(csv.DictReader(f))
    with args.native_units.open(newline="", encoding="utf-8") as f:
        native_rows = list(csv.DictReader(f))
    with args.wgsrpd_support.open(newline="", encoding="utf-8") as f:
        support_rows = list(csv.DictReader(f))

    if len(candidates) != 642 or len({r["species"] for r in candidates}) != 642:
        raise RuntimeError("frozen 642 species identity drift")
    verify_host_pair_identity(candidates, pair_rows)

    assets = {}
    for path in args.current_asset:
        name = path.name
        matched = [v for v in VARIABLE_ORDER if name == f"CHELSA_TraCE21K_{v}_20_V1.0.tif"]
        if len(matched) != 1 or matched[0] in assets:
            raise RuntimeError(f"unfrozen or duplicate CHELSA current asset: {name}")
        assets[matched[0]] = path
    if set(assets) != set(VARIABLE_ORDER):
        raise RuntimeError("CHELSA 0-BP required variables not complete")
    coords = [(float(row["longitude"]), float(row["latitude"])) for row in support_rows]
    values = sample_matrix(coords, [assets[v] for v in VARIABLE_ORDER])
    finite = np.all(np.isfinite(values), axis=1)
    if sum(finite) < 100:
        raise RuntimeError("frozen reference cloud has fewer than 100 finite points")
    # No response, no beta, and no inferred paleo occurrence are read.
    counts = census_finite_host_support(candidates, native_rows, support_rows, finite)
    write_receipt(args.candidates, args.host_pairs, args.native_units,
                  args.wgsrpd_support, [assets[v] for v in VARIABLE_ORDER],
                  counts, sum(finite), args.rule, args.output_receipt, args.output_species)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
