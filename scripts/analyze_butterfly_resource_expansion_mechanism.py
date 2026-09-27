#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from collections import defaultdict
from pathlib import Path

import numpy as np


EXPECTED_DESCRIPTOR_SHA256 = (
    "894f48dbca1760fc4fa75bfee8f663540ab4b9380f8b9daf2bc09440e8bb0cdc"
)


def sha256_path(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


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


def host_band(value: float) -> str:
    if value == 1:
        return "1_family"
    if value == 2:
        return "2_families"
    if 3 <= value <= 5:
        return "3_to_5_families"
    if value >= 6:
        return "6plus_families"
    raise ValueError(f"invalid host-family count: {value}")


def fractional_host_contributions(
    host_added_units: dict[str, set[str]],
    butterfly_added_union: set[str],
) -> dict[str, object]:
    if not butterfly_added_union:
        return {
            "introduced_added_units": 0,
            "contributing_host_species": 0,
            "maximum_single_host_fractional_share": None,
            "top_two_host_fractional_share": None,
            "effective_contributor_number": None,
            "fraction_added_units_with_multiple_contributing_hosts": None,
            "fractional_credits": {},
        }

    credit: dict[str, float] = defaultdict(float)
    multiple = 0
    for unit in sorted(butterfly_added_union):
        contributors = sorted(
            host
            for host, added in host_added_units.items()
            if unit in added
        )
        if not contributors:
            raise RuntimeError(
                f"added butterfly unit {unit} has no contributing host species"
            )
        if len(contributors) > 1:
            multiple += 1
        share = 1.0 / len(contributors)
        for host in contributors:
            credit[host] += share

    total = float(len(butterfly_added_union))
    if abs(sum(credit.values()) - total) > 1e-9:
        raise RuntimeError("fractional host credits do not sum to added union")

    proportions = sorted(
        (value / total for value in credit.values()),
        reverse=True,
    )
    max_share = proportions[0]
    top_two = sum(proportions[:2])
    effective = 1.0 / sum(value * value for value in proportions)

    return {
        "introduced_added_units": int(total),
        "contributing_host_species": len(proportions),
        "maximum_single_host_fractional_share": float(max_share),
        "top_two_host_fractional_share": float(top_two),
        "effective_contributor_number": float(effective),
        "fraction_added_units_with_multiple_contributing_hosts": (
            multiple / total
        ),
        "fractional_credits": {
            host: float(value)
            for host, value in sorted(credit.items())
        },
    }


def load_descriptors(path: Path) -> dict[str, dict[str, object]]:
    if sha256_path(path) != EXPECTED_DESCRIPTOR_SHA256:
        raise RuntimeError("descriptor SHA drift")
    rows = {}
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        required = {
            "species",
            "host_family_count",
            "host_wgsrpd3_unit_count",
            "resolved_host_species",
        }
        if not required <= set(reader.fieldnames or ()):
            raise RuntimeError("descriptor schema drift")
        for row in reader:
            name = str(row["species"]).strip()
            rows[name] = {
                "species": name,
                "host_family_count": float(row["host_family_count"]),
                "native_resource_units_descriptor": int(
                    row["host_wgsrpd3_unit_count"]
                ),
                "resolved_host_species": int(row["resolved_host_species"]),
            }
    if len(rows) != 339:
        raise RuntimeError("expected exact 339-species descriptor table")
    return rows


def load_pairs(path: Path) -> dict[str, set[str]]:
    out: dict[str, set[str]] = defaultdict(set)
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        required = {"insect_species", "accepted_plant_name_id"}
        if not required <= set(reader.fieldnames or ()):
            raise RuntimeError("interaction sidecar schema drift")
        for row in reader:
            insect = str(row["insect_species"]).strip()
            host = str(row["accepted_plant_name_id"]).strip()
            if insect and host:
                out[insect].add(host)
    return out


def load_units(path: Path) -> dict[str, set[str]]:
    out: dict[str, set[str]] = defaultdict(set)
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        required = {"accepted_plant_name_id", "area_code_l3"}
        if not required <= set(reader.fieldnames or ()):
            raise RuntimeError("distribution sidecar schema drift")
        for row in reader:
            host = str(row["accepted_plant_name_id"]).strip()
            unit = str(row["area_code_l3"]).strip()
            if host and unit:
                out[host].add(unit)
    return out


def percentile(values, q: float) -> float | None:
    values = list(values)
    if not values:
        return None
    return float(np.quantile(np.asarray(values, dtype=float), q))


def summarize(rows: list[dict[str, object]]) -> dict[str, object]:
    expanded = [row for row in rows if int(row["introduced_added_units"]) > 0]
    host_family = [float(row["host_family_count"]) for row in expanded]
    max_share = [
        float(row["maximum_single_host_fractional_share"])
        for row in expanded
    ]
    effective = [
        float(row["effective_contributor_number"])
        for row in expanded
    ]
    resolved = [int(row["resolved_host_species"]) for row in expanded]

    bands = {}
    for label in ("1_family", "2_families", "3_to_5_families", "6plus_families"):
        subset = [
            row
            for row in expanded
            if row["host_breadth_stratum"] == label
        ]
        if not subset:
            continue
        bands[label] = {
            "expanded_species": len(subset),
            "median_maximum_single_host_fractional_share": percentile(
                [row["maximum_single_host_fractional_share"] for row in subset],
                0.5,
            ),
            "median_top_two_host_fractional_share": percentile(
                [row["top_two_host_fractional_share"] for row in subset],
                0.5,
            ),
            "median_effective_contributor_number": percentile(
                [row["effective_contributor_number"] for row in subset],
                0.5,
            ),
            "median_contributing_host_species": percentile(
                [row["contributing_host_species"] for row in subset],
                0.5,
            ),
            "median_fraction_multi_host_added_units": percentile(
                [
                    row["fraction_added_units_with_multiple_contributing_hosts"]
                    for row in subset
                ],
                0.5,
            ),
        }

    return {
        "species": len(rows),
        "expanded_species": len(expanded),
        "maximum_single_host_fractional_share": {
            "q25": percentile(max_share, 0.25),
            "median": percentile(max_share, 0.5),
            "q75": percentile(max_share, 0.75),
        },
        "top_two_host_fractional_share": {
            "median": percentile(
                [row["top_two_host_fractional_share"] for row in expanded],
                0.5,
            )
        },
        "effective_contributor_number": {
            "q25": percentile(effective, 0.25),
            "median": percentile(effective, 0.5),
            "q75": percentile(effective, 0.75),
        },
        "fraction_added_units_with_multiple_contributing_hosts": {
            "median": percentile(
                [
                    row["fraction_added_units_with_multiple_contributing_hosts"]
                    for row in expanded
                ],
                0.5,
            )
        },
        "spearman": {
            "host_family_count_vs_maximum_single_host_fractional_share": spearman(
                host_family, max_share
            ),
            "host_family_count_vs_effective_contributor_number": spearman(
                host_family, effective
            ),
            "resolved_host_species_vs_effective_contributor_number": spearman(
                resolved, effective
            ),
        },
        "host_breadth_strata": bands,
    }


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--protocol-json", type=Path, required=True)
    ap.add_argument("--descriptors-csv", type=Path, required=True)
    ap.add_argument("--insect-host-csv", type=Path, required=True)
    ap.add_argument("--native-distribution-csv", type=Path, required=True)
    ap.add_argument("--contemporary-distribution-csv", type=Path, required=True)
    ap.add_argument("--output-csv", type=Path, required=True)
    ap.add_argument("--output-json", type=Path, required=True)
    args = ap.parse_args()

    protocol = json.loads(args.protocol_json.read_text(encoding="utf-8"))
    if protocol.get("schema") != (
        "ttf_butterfly_resource_expansion_mechanism_protocol_v0.1"
    ):
        raise RuntimeError("unexpected mechanism protocol schema")
    if protocol.get("status") != (
        "FROZEN_EXPLORATORY_AFTER_AGGREGATE_RESULT_BEFORE_HOST_CONTRIBUTION_DECOMPOSITION"
    ):
        raise RuntimeError("mechanism protocol is not frozen")

    descriptors = load_descriptors(args.descriptors_csv)
    pairs = load_pairs(args.insect_host_csv)
    native_by_host = load_units(args.native_distribution_csv)
    contemporary_by_host = load_units(args.contemporary_distribution_csv)

    rows = []
    host_credit_payload = {}
    for name, descriptor in descriptors.items():
        host_families = float(descriptor["host_family_count"])
        if (
            host_families <= 0
            or int(descriptor["native_resource_units_descriptor"]) <= 0
        ):
            continue

        hosts = sorted(pairs.get(name, set()))
        if not hosts:
            raise RuntimeError(f"resource-eligible species lacks resolved hosts: {name}")

        native_union: set[str] = set()
        contemporary_union: set[str] = set()
        added_by_host: dict[str, set[str]] = {}
        for host in hosts:
            native_units = set(native_by_host.get(host, set()))
            contemporary_units = set(contemporary_by_host.get(host, set()))
            if not native_units.issubset(contemporary_units):
                raise RuntimeError(f"native host units not subset contemporary: {host}")
            native_union.update(native_units)
            contemporary_union.update(contemporary_units)
            added_by_host[host] = contemporary_units - native_units

        if len(native_union) != int(descriptor["native_resource_units_descriptor"]):
            raise RuntimeError(
                f"native union drift for {name}: "
                f"{len(native_union)} vs {descriptor['native_resource_units_descriptor']}"
            )

        added_union = contemporary_union - native_union
        contribution = fractional_host_contributions(
            added_by_host,
            added_union,
        )
        row = {
            "species": name,
            "host_family_count": host_families,
            "host_breadth_stratum": host_band(host_families),
            "resolved_host_species": int(descriptor["resolved_host_species"]),
            "native_resource_units": len(native_union),
            "contemporary_resource_units": len(contemporary_union),
            **{
                key: value
                for key, value in contribution.items()
                if key != "fractional_credits"
            },
            "host_taxonomy_lower_bound_adequate": (
                int(descriptor["resolved_host_species"]) >= host_families
            ),
        }
        rows.append(row)
        if added_union:
            host_credit_payload[name] = contribution["fractional_credits"]

    rows = sorted(rows, key=lambda row: str(row["species"]))
    if len(rows) != 239:
        raise RuntimeError(f"expected 239 resource-eligible species; got {len(rows)}")

    adequate = [
        row
        for row in rows
        if bool(row["host_taxonomy_lower_bound_adequate"])
    ]
    payload = {
        "schema": "ttf_butterfly_resource_expansion_mechanism_result_v0.1",
        "status": "EXPLORATORY_HOST_CONTRIBUTION_DECOMPOSITION_COMPLETE",
        "protocol": str(args.protocol_json),
        "full_resource_eligible_panel": summarize(rows),
        "host_taxonomy_lower_bound_adequate_subset": summarize(adequate),
        "host_fractional_credits_for_expanded_species": host_credit_payload,
        "claim_boundary": protocol["claim_boundary"],
    }

    args.output_csv.parent.mkdir(parents=True, exist_ok=True)
    fields = list(rows[0].keys())
    with args.output_csv.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)
    args.output_json.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    print(
        json.dumps(
            {
                "schema": payload["schema"],
                "status": payload["status"],
                "full_resource_eligible_panel": payload[
                    "full_resource_eligible_panel"
                ],
                "host_taxonomy_lower_bound_adequate_subset": payload[
                    "host_taxonomy_lower_bound_adequate_subset"
                ],
            },
            indent=2,
            sort_keys=True,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
