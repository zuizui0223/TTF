"""Exact source-only sensitivity of Ellis 2007 quartz identifiability to restoring original arms.

Uses only the frozen 24-vs-18 topology receipt. Exhausts all 2^6 subsets
of the six original within-basin cross-quartz directed transplant arms.
No fitness, genetics, GPS, or historical-climate data are read.
"""
from __future__ import annotations
import argparse
import hashlib
import itertools
import json
from collections import Counter
from pathlib import Path

from ellis2007_design_rank_v01 import RECEIPT, SHA, analyze, rank, category, CLASSES

BASINS = ("Ariz", "Sout", "QK")
QUARTZ = ("dense", "sparse")


def row(source, garden):
    return [
        1,
        int(source[1] == "sparse"),
        int(garden[1] == "sparse"),
        int(source[0] == "Sout"),
        int(source[0] == "QK"),
        int(garden[0] == "Sout"),
        int(garden[0] == "QK"),
    ]


def components(nodes, edges):
    adjacency = {x: set() for x in nodes}
    for a, b in edges:
        adjacency[a].add(b)
        adjacency[b].add(a)
    remaining = set(nodes)
    groups = []
    while remaining:
        pending = [min(remaining)]
        seen = set()
        while pending:
            u = pending.pop()
            if u in seen:
                continue
            seen.add(u)
            pending.extend(adjacency[u] - seen)
        groups.append(seen)
        remaining -= seen
    return groups


def evaluate(data):
    analyze(data)
    nodes = [(basin, quartz) for basin in BASINS for quartz in QUARTZ]
    secondary = [(s, g) for s in nodes for g in nodes
                 if category(s, g) in (CLASSES[0], CLASSES[2])]
    missing = [(s, g) for s in nodes for g in nodes if category(s, g) == CLASSES[1]]
    if len(secondary) != 18 or len(missing) != 6:
        raise AssertionError("original source-only topology changed")
    if rank([row(s, g) for s, g in secondary]) != 6:
        raise AssertionError("frozen secondary rank changed")

    results = {}
    for count in range(7):
        designs = []
        for added in itertools.combinations(missing, count):
            edges = secondary + list(added)
            designs.append({
                "rank": rank([row(s, g) for s, g in edges]),
                "components": len(components(nodes, edges)),
                "basin_coverage": len({s[0] for s, g in added}),
            })
        results[str(count)] = {
            "subsets": len(designs),
            "full_rank_subsets": sum(x["rank"] == 7 for x in designs),
            "connected_subsets": sum(x["components"] == 1 for x in designs),
            "min_basin_coverage": min(x["basin_coverage"] for x in designs),
            "max_basin_coverage": max(x["basin_coverage"] for x in designs),
        }
    if [results[str(k)]["subsets"] for k in range(7)] != [1, 6, 15, 20, 15, 6, 1]:
        raise AssertionError("exhaustive six-arm subset census drift")
    if results["0"]["full_rank_subsets"] != 0 or any(
        results[str(k)]["full_rank_subsets"] != results[str(k)]["subsets"] for k in range(1, 7)
    ):
        raise AssertionError("restoration rank prediction failed")
    if any(results[str(k)]["connected_subsets"] != results[str(k)]["subsets"] for k in range(1, 7)):
        raise AssertionError("restoration connectivity prediction failed")

    return {
        "schema": "ttf_Ellis2007_source_only_missing_arm_restoration_v01",
        "source_receipt_sha256": SHA,
        "status": "DESIGN_RESTORATION_IDENTIFIABILITY_ONLY_NOT_EMPIRICAL_RECOVERY",
        "secondary_directed_arms": 18,
        "original_missing_within_basin_cross_quartz_arms": 6,
        "subset_count": 64,
        "baseline_rank": 6,
        "baseline_components": 2,
        "minimum_added_original_arms_for_full_rank": 1,
        "minimum_added_original_arms_for_connected_graph": 1,
        "by_number_of_restored_arms": results,
        "all_six_arms_present_in_original_experimental_design": True,
        "actual_primary_response_values_recovered": False,
        "restored_one_arm_sufficient_for_precision_or_causal_identification": False,
        "original_source_gps_verified": False,
        "fitness_genetics_LGM_or_adaptation_response_opened": False,
        "frozen_TTF_endpoints_unchanged": True,
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    raw = RECEIPT.read_bytes()
    if hashlib.sha256(raw).hexdigest() != SHA:
        raise ValueError("frozen source receipt SHA mismatch")
    result = evaluate(json.loads(raw))
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, sort_keys=True, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(result, sort_keys=True))


if __name__ == "__main__":
    main()
