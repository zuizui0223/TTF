#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import numpy as np


STRATA = ("1_family", "2_families", "3_to_5_families", "6plus_families")


def average_ranks(values) -> np.ndarray:
    x = np.asarray(values, dtype=float)
    order = np.argsort(x, kind="mergesort")
    ranks = np.empty(len(x), dtype=float)
    start = 0
    while start < len(order):
        stop = start + 1
        while stop < len(order) and x[order[stop]] == x[order[start]]:
            stop += 1
        ranks[order[start:stop]] = 0.5 * ((start + 1) + stop)
        start = stop
    return ranks


def spearman(a, b) -> float | None:
    if len(a) < 3 or len(a) != len(b):
        return None
    x = average_ranks(a)
    y = average_ranks(b)
    if np.std(x) <= np.sqrt(np.finfo(float).eps):
        return None
    if np.std(y) <= np.sqrt(np.finfo(float).eps):
        return None
    return float(np.corrcoef(x, y)[0, 1])


def richness_tertiles(rows: list[dict[str, object]]) -> dict[str, dict[str, float | int]]:
    ordered = sorted(
        rows,
        key=lambda row: (
            int(row["resolved_host_species"]),
            str(row["species"]),
        ),
    )
    n = len(ordered)
    if n < 3:
        return {}
    groups = {
        "low": ordered[: n // 3],
        "mid": ordered[n // 3 : 2 * n // 3],
        "high": ordered[2 * n // 3 :],
    }
    out = {}
    for label, subset in groups.items():
        if not subset:
            continue
        out[label] = {
            "species": len(subset),
            "resolved_host_species_min": min(
                int(row["resolved_host_species"]) for row in subset
            ),
            "resolved_host_species_median": float(
                np.median([int(row["resolved_host_species"]) for row in subset])
            ),
            "resolved_host_species_max": max(
                int(row["resolved_host_species"]) for row in subset
            ),
            "median_effective_contributor_number": float(
                np.median(
                    [float(row["effective_contributor_number"]) for row in subset]
                )
            ),
            "median_maximum_single_host_fractional_share": float(
                np.median(
                    [
                        float(row["maximum_single_host_fractional_share"])
                        for row in subset
                    ]
                )
            ),
        }
    return out


def summarize_stratum(rows: list[dict[str, object]]) -> dict[str, object]:
    richness = [int(row["resolved_host_species"]) for row in rows]
    effective = [float(row["effective_contributor_number"]) for row in rows]
    max_share = [
        float(row["maximum_single_host_fractional_share"]) for row in rows
    ]
    added = [int(row["introduced_added_units"]) for row in rows]
    return {
        "species": len(rows),
        "resolved_host_species": {
            "minimum": min(richness),
            "median": float(np.median(richness)),
            "maximum": max(richness),
        },
        "spearman": {
            "resolved_host_species_vs_effective_contributor_number": spearman(
                richness, effective
            ),
            "resolved_host_species_vs_maximum_single_host_fractional_share": spearman(
                richness, max_share
            ),
            "resolved_host_species_vs_introduced_added_units": spearman(
                richness, added
            ),
        },
        "resolved_host_species_tertiles": richness_tertiles(rows),
    }


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--protocol-json", type=Path, required=True)
    ap.add_argument("--species-metrics-csv", type=Path, required=True)
    ap.add_argument("--output-json", type=Path, required=True)
    args = ap.parse_args()

    protocol = json.loads(args.protocol_json.read_text(encoding="utf-8"))
    if protocol.get("schema") != "ttf_butterfly_host_specialization_hierarchy_protocol_v0.1":
        raise RuntimeError("unexpected hierarchy protocol schema")
    if protocol.get("status") != (
        "FROZEN_EXPLORATORY_AFTER_RESOURCE_ARCHITECTURE_RESULT_BEFORE_WITHIN_STRATUM_ASSOCIATIONS"
    ):
        raise RuntimeError("hierarchy protocol is not frozen")

    rows = []
    with args.species_metrics_csv.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        required = {
            "species",
            "host_breadth_stratum",
            "resolved_host_species",
            "introduced_added_units",
            "effective_contributor_number",
            "maximum_single_host_fractional_share",
            "host_taxonomy_lower_bound_adequate",
        }
        if not required <= set(reader.fieldnames or ()):
            raise RuntimeError("mechanism metric schema drift")
        for row in reader:
            adequate = str(row["host_taxonomy_lower_bound_adequate"]).strip().lower()
            if adequate not in {"true", "1"}:
                continue
            if int(row["introduced_added_units"]) <= 0:
                continue
            rows.append(
                {
                    "species": str(row["species"]).strip(),
                    "host_breadth_stratum": str(
                        row["host_breadth_stratum"]
                    ).strip(),
                    "resolved_host_species": int(row["resolved_host_species"]),
                    "introduced_added_units": int(row["introduced_added_units"]),
                    "effective_contributor_number": float(
                        row["effective_contributor_number"]
                    ),
                    "maximum_single_host_fractional_share": float(
                        row["maximum_single_host_fractional_share"]
                    ),
                }
            )

    strata = {}
    for label in STRATA:
        subset = [
            row for row in rows if row["host_breadth_stratum"] == label
        ]
        if len(subset) < 3:
            raise RuntimeError(f"stratum {label} has too few rows: {len(subset)}")
        strata[label] = summarize_stratum(subset)

    payload = {
        "schema": "ttf_butterfly_host_specialization_hierarchy_result_v0.1",
        "status": "EXPLORATORY_WITHIN_STRATUM_ASSOCIATIONS_COMPLETE",
        "eligible_expanded_species": len(rows),
        "primary_focus": {
            "stratum": "1_family",
            **strata["1_family"],
        },
        "strata": strata,
        "interpretation": {
            "family_breadth_held_constant_in_primary_focus": True,
            "question": (
                "Whether species-level host richness predicts portfolio-style "
                "anthropogenic resource opportunity even among butterflies that "
                "all use exactly one host family."
            ),
        },
        "claim_boundary": protocol["claim_boundary"],
    }
    args.output_json.parent.mkdir(parents=True, exist_ok=True)
    args.output_json.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    print(json.dumps(payload, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
