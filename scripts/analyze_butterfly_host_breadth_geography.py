#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path

import numpy as np
def load_rows(path: Path) -> list[dict[str, object]]:
    rows = []
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        required = {
            "species",
            "host_family_count",
            "host_wgsrpd3_unit_count",
            "resolved_host_species",
            "hosts_with_primary_native_units",
        }
        if not required <= set(reader.fieldnames or ()):
            raise RuntimeError("resource descriptor schema drift")
        for row in reader:
            rows.append(
                {
                    "species": str(row["species"]),
                    "host_family_count": float(row["host_family_count"]),
                    "host_wgsrpd3_unit_count": int(row["host_wgsrpd3_unit_count"]),
                    "resolved_host_species": int(row["resolved_host_species"]),
                    "hosts_with_primary_native_units": int(
                        row["hosts_with_primary_native_units"]
                    ),
                }
            )
    return rows


def band(value: float) -> str:
    if value == 1:
        return "1_family"
    if value == 2:
        return "2_families"
    if value <= 5:
        return "3_to_5_families"
    return "6plus_families"


def quantile(values: list[float], q: float) -> float:
    return float(np.quantile(np.asarray(values, dtype=float), q))


def average_ranks(values: np.ndarray) -> np.ndarray:
    values = np.asarray(values, dtype=float)
    order = np.argsort(values, kind="mergesort")
    ranks = np.empty(len(values), dtype=float)
    start = 0
    while start < len(order):
        stop = start + 1
        while stop < len(order) and values[order[stop]] == values[order[start]]:
            stop += 1
        average = 0.5 * ((start + 1) + stop)
        ranks[order[start:stop]] = average
        start = stop
    return ranks


def spearman(values_a: np.ndarray, values_b: np.ndarray) -> float:
    a = average_ranks(np.asarray(values_a, dtype=float))
    b = average_ranks(np.asarray(values_b, dtype=float))
    if np.std(a) == 0 or np.std(b) == 0:
        return float("nan")
    return float(np.corrcoef(a, b)[0, 1])


def analyze(rows: list[dict[str, object]]) -> dict:
    eligible = [
        row
        for row in rows
        if float(row["host_family_count"]) > 0
        and int(row["host_wgsrpd3_unit_count"]) > 0
    ]
    family = np.asarray([float(r["host_family_count"]) for r in eligible], float)
    units = np.asarray([int(r["host_wgsrpd3_unit_count"]) for r in eligible], float)
    resolved = np.asarray([int(r["resolved_host_species"]) for r in eligible], float)

    rho_family_units = spearman(family, units)
    rho_resolved_units = spearman(resolved, units)
    rho_family_resolved = spearman(family, resolved)

    lower_bound_adequate = [
        r
        for r in eligible
        if int(r["resolved_host_species"]) >= float(r["host_family_count"])
    ]
    adequate_rho = spearman(
        np.asarray(
            [float(r["host_family_count"]) for r in lower_bound_adequate],
            dtype=float,
        ),
        np.asarray(
            [int(r["host_wgsrpd3_unit_count"]) for r in lower_bound_adequate],
            dtype=float,
        ),
    )

    y = np.log1p(units)
    x = np.column_stack(
        [
            np.ones(len(eligible)),
            np.log1p(family),
            np.log1p(resolved),
        ]
    )
    coef, *_ = np.linalg.lstsq(x, y, rcond=None)
    fitted = x @ coef
    ss_total = float(np.sum((y - np.mean(y)) ** 2))
    ss_resid = float(np.sum((y - fitted) ** 2))
    r2 = None if ss_total <= 0 else 1.0 - ss_resid / ss_total

    bands = {}
    for label in ("1_family", "2_families", "3_to_5_families", "6plus_families"):
        subset = [
            r
            for r in eligible
            if band(float(r["host_family_count"])) == label
        ]
        vals = [int(r["host_wgsrpd3_unit_count"]) for r in subset]
        bands[label] = {
            "species": len(subset),
            "median_host_wgsrpd3_units": quantile(vals, 0.5),
            "q25_host_wgsrpd3_units": quantile(vals, 0.25),
            "q75_host_wgsrpd3_units": quantile(vals, 0.75),
        }

    specialist_wide = sorted(
        (
            r
            for r in eligible
            if float(r["host_family_count"]) == 1
        ),
        key=lambda r: (-int(r["host_wgsrpd3_unit_count"]), str(r["species"])),
    )[:10]
    generalist_narrow = sorted(
        (
            r
            for r in eligible
            if float(r["host_family_count"]) >= 6
        ),
        key=lambda r: (int(r["host_wgsrpd3_unit_count"]), str(r["species"])),
    )[:10]

    return {
        "schema": "ttf_butterfly_host_breadth_geography_exploratory_v0.1",
        "status": "EXPLORATORY_DESCRIPTIVE_RESULT",
        "species": {
            "S1_total": len(rows),
            "resource_eligible": len(eligible),
            "lower_bound_host_taxonomy_adequate": len(lower_bound_adequate),
            "lower_bound_definition": (
                "HOSTS-WCVP resolved host-species count >= LepTraits "
                "NumberOfHostplantFamilies; necessary but not sufficient for "
                "host-record completeness"
            ),
        },
        "rank_correlations": {
            "host_family_count_vs_native_host_wgsrpd3_units": rho_family_units,
            "resolved_host_species_vs_native_host_wgsrpd3_units": rho_resolved_units,
            "host_family_count_vs_resolved_host_species": rho_family_resolved,
            "host_family_count_vs_native_host_wgsrpd3_units_lower_bound_adequate_only": adequate_rho,
        },
        "log_linear_descriptive_model": {
            "response": "log1p(native host-resource WGSRPD3 unit count)",
            "predictors": [
                "intercept",
                "log1p(LepTraits host-family count)",
                "log1p(HOSTS-WCVP resolved host-species count)",
            ],
            "coefficients": [float(v) for v in coef],
            "r_squared": None if r2 is None else float(r2),
            "inference": "descriptive only; no confirmatory p-value",
        },
        "host_family_bands": bands,
        "counterexamples": {
            "ten_widest_resource_footprints_among_one_family_species": specialist_wide,
            "ten_narrowest_resource_footprints_among_sixplus_family_species": generalist_narrow,
        },
        "interpretation": (
            "Taxonomic host-family breadth and geographic host-resource breadth "
            "are related but are not interchangeable. Much of the geographic "
            "breadth signal tracks the number of resolved host species rather "
            "than host-family count itself, and strong specialist/generalist "
            "counterexamples occur."
        ),
        "claim_boundary": {
            "exploratory": True,
            "host_database_completeness_is_not_proven": True,
            "native_host_geography_only": True,
            "butterfly_occurrences_used": False,
            "genetic_response_used": False,
        },
    }


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--descriptors-csv", type=Path, required=True)
    ap.add_argument("--output-json", type=Path, required=True)
    args = ap.parse_args()
    result = analyze(load_rows(args.descriptors_csv))
    args.output_json.parent.mkdir(parents=True, exist_ok=True)
    args.output_json.write_text(
        json.dumps(result, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    print(json.dumps(result, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
