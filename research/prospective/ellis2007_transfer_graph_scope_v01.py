"""Describe which Ellis 2007 transplant comparison classes survive curation.

Source-only and descriptive: the only input is an immutable prior aggregate
receipt. No raw fitness, phenotype, AFLP, climate, GPS or specimen IDs.
"""
from __future__ import annotations
import argparse
import hashlib
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
RECEIPT = HERE / "source_registry/ellis2007_24_to_18_origin_garden_arm_audit_v02.json"
SHA256 = "28cea27d59affd846c9ee490f20daedd36220e6e75041f429663a208dcd239ce"
BASINS = ("Ariz", "Sout", "QK")
QUARTZ = ("dense", "sparse")
ARM_CLASSES = ("local_same_source_and_garden", "same_basin_other_quartz",
               "other_basin_same_quartz", "other_basin_other_quartz_not_in_paper_design")


def arm_class(source, garden):
    if source == garden:
        return ARM_CLASSES[0]
    if source[0] == garden[0]:
        return ARM_CLASSES[1]
    if source[1] == garden[1]:
        return ARM_CLASSES[2]
    return ARM_CLASSES[3]


def connected_components(nodes, edges):
    links = {n: set() for n in nodes}
    for source, garden in edges:
        links[source].add(garden)
        links[garden].add(source)
    remaining = set(nodes)
    answer = []
    while remaining:
        frontier, seen = [min(remaining)], set()
        while frontier:
            node = frontier.pop()
            if node in seen:
                continue
            seen.add(node)
            frontier.extend(links[node] - seen)
        remaining -= seen
        answer.append(seen)
    return answer


def analyze(receipt):
    if receipt.get("schema") != "ttf_Ellis2007_primary_24_source_garden_vs_Hargreaves_18_source_only_edges_diagnosis_v02":
        raise ValueError("wrong frozen receipt schema")
    if receipt.get("status") != "PRIMARY_ORIGINAL_24_SOURCE_GARDEN_DESIGN_VS_SECONDARY_EXTRACT_18_WITH_SOURCE_ARM_MISSINGNESS":
        raise ValueError("wrong source-only frozen status")
    for field in ("no_fitness_growth_survival_reproduction_or_AFLP_response_rows_read",
                  "no_LGM_21ka_0BP_climate_effect_estimated",
                  "old_Bontrager_44of50_terminal_stop_unchanged"):
        if receipt.get(field) is not True:
            raise ValueError("frozen protection changed: " + field)
    rows = receipt.get("published_origin_site_pairs_by_class")
    if not isinstance(rows, list) or len(rows) != 4:
        raise ValueError("invalid frozen rows")
    by_class = {r["class"]: r for r in rows}
    if len(by_class) != 4 or set(by_class) != set(ARM_CLASSES):
        raise ValueError("duplicate or unknown arm class")
    expected, observed = {}, {}
    for label in ARM_CLASSES:
        row = by_class[label]
        p, o, missing = (row["paper_expected"],
                         row["original_Hargreaves_extract_observed"],
                         row["not_represented_in_secondary_original"])
        if any(type(v) is not int or v < 0 for v in (p, o, missing)) or p != o + missing:
            raise ValueError("invalid class census")
        expected[label], observed[label] = p, o
    if tuple(expected.values()) != (6, 6, 12, 0):
        raise ValueError("primary design drift")
    if (sum(expected.values()) != receipt["original_primary_expected_directed_source_garden_pairs"]
            or sum(observed.values()) != receipt["observed_Hargreaves_distinct_directed_source_garden_pairs"]):
        raise ValueError("class sum and total disagree")

    nodes = [(b, h) for b in BASINS for h in QUARTZ]
    class_edges = {k: set() for k in ARM_CLASSES}
    for source in nodes:
        for garden in nodes:
            class_edges[arm_class(source, garden)].add((source, garden))
    assert tuple(map(lambda c: len(class_edges[c]), ARM_CLASSES)) == (6, 6, 12, 12)
    # Aggregate data cannot identify which individual edges survived if only
    # PART of a class survived. Stop instead of inventing raw experiment IDs.
    if any(observed[k] not in (0, len(class_edges[k])) for k in ARM_CLASSES):
        raise ValueError("partial class: individual edges unidentifiable")
    p_edges = set().union(*(class_edges[k] for k in ARM_CLASSES if expected[k]))
    s_edges = set().union(*(class_edges[k] for k in ARM_CLASSES if observed[k]))
    p_parts = connected_components(nodes, p_edges)
    s_parts = connected_components(nodes, s_edges)
    if len(p_parts) != 1 or sorted(map(len, s_parts)) != [3, 3]:
        raise AssertionError("unexpected source-only comparison connectivity")
    if not all(len({h for _, h in part}) == 1 for part in s_parts):
        raise AssertionError("secondary comparison graph does not split by quartz")
    return {
        "schema": "ttf_Ellis2007_source_only_transfer_graph_scope_v01",
        "source_receipt_sha256": SHA256,
        "status": "DESCRIPTIVE_COMPARISON_SUPPORT_NOT_BIOLOGICAL_OR_LGM_EFFECT",
        "original_directed_arms": len(p_edges),
        "secondary_directed_arms": len(s_edges),
        "original_connected_components": len(p_parts),
        "secondary_connected_components": len(s_parts),
        "secondary_component_sizes": sorted(map(len, s_parts)),
        "secondary_components_homogeneous_by_quartz": True,
        "original_within_basin_cross_quartz_arms": expected[ARM_CLASSES[1]],
        "secondary_within_basin_cross_quartz_arms": observed[ARM_CLASSES[1]],
        "original_observed_arm_classes": sum(v > 0 for v in expected.values()),
        "secondary_observed_arm_classes": sum(v > 0 for v in observed.values()),
        "secondary_local_vs_other_basin_same_quartz_comparison_present": True,
        "secondary_within_basin_cross_quartz_comparison_present": False,
        "genetic_or_fitness_or_LGM_response_opened": False,
        "former_frozen_endpoints_changed": False,
        "source_garden_IDs_or_coordinates_exported": False,
        "graph_edges_are_transplants_not_genetic_migration": True,
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    raw = RECEIPT.read_bytes()
    if hashlib.sha256(raw).hexdigest() != SHA256:
        raise ValueError("frozen original receipt SHA mismatch")
    data = analyze(json.loads(raw))
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(data, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(data, sort_keys=True))


if __name__ == "__main__":
    main()
