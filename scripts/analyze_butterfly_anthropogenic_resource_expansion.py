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

from ttf.lepidoptera_host_resource import build_insect_host_footprints


EXPECTED_DESCRIPTOR_SHA256 = (
    "894f48dbca1760fc4fa75bfee8f663540ab4b9380f8b9daf2bc09440e8bb0cdc"
)
EXPECTED_NATIVE_FOOTPRINT_SHA256 = (
    "546c38970789a153735d6e54cc66bdd2e8987a22b93b097aa2062e8c7f84baef"
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
    if len(a) < 3 or len(b) != len(a):
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
    raise ValueError(f"invalid host family count: {value}")


def load_descriptors(path: Path) -> dict[str, dict[str, object]]:
    if sha256_path(path) != EXPECTED_DESCRIPTOR_SHA256:
        raise RuntimeError("descriptor SHA-256 drift")
    rows = {}
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        required = {
            "species",
            "host_family_count",
            "host_wgsrpd3_unit_count",
            "resolved_host_species",
            "hosts_with_primary_native_units",
            "wing_size_proxy",
            "voltinism",
        }
        if not required <= set(reader.fieldnames or ()):
            raise RuntimeError("descriptor schema drift")
        for row in reader:
            name = str(row["species"]).strip()
            rows[name] = {
                "species": name,
                "host_family_count": float(row["host_family_count"]),
                "native_descriptor_units": int(row["host_wgsrpd3_unit_count"]),
                "resolved_host_species": int(row["resolved_host_species"]),
                "hosts_with_primary_native_units": int(
                    row["hosts_with_primary_native_units"]
                ),
                "wing_size_proxy": float(row["wing_size_proxy"]),
                "voltinism": str(row["voltinism"]).strip(),
            }
    if len(rows) != 339:
        raise RuntimeError(f"expected 339 S1 descriptor species; got {len(rows)}")
    return rows


def load_native_footprints(path: Path) -> dict[str, frozenset[str]]:
    if sha256_path(path) != EXPECTED_NATIVE_FOOTPRINT_SHA256:
        raise RuntimeError("native footprint SHA-256 drift")
    payload = json.loads(path.read_text(encoding="utf-8"))
    if payload.get("schema") != "ttf_butterfly_resource_envelope_s1_footprints_v0.1":
        raise RuntimeError("unexpected native footprint schema")
    return {
        str(name): frozenset(map(str, units))
        for name, units in payload.get("species", {}).items()
    }


def load_contemporary_sidecars(
    interaction_path: Path,
    distribution_path: Path,
):
    pairs: list[tuple[str, str]] = []
    with interaction_path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        required = {"insect_species", "accepted_plant_name_id"}
        if not required <= set(reader.fieldnames or ()):
            raise RuntimeError("insect-host sidecar schema drift")
        for row in reader:
            insect = str(row["insect_species"]).strip()
            host = str(row["accepted_plant_name_id"]).strip()
            if insect and host:
                pairs.append((insect, host))

    units: dict[str, set[str]] = defaultdict(set)
    with distribution_path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        required = {"accepted_plant_name_id", "area_code_l3"}
        if not required <= set(reader.fieldnames or ()):
            raise RuntimeError("contemporary distribution sidecar schema drift")
        for row in reader:
            host = str(row["accepted_plant_name_id"]).strip()
            unit = str(row["area_code_l3"]).strip()
            if host and unit:
                units[host].add(unit)

    return build_insect_host_footprints(pairs, units)


def descriptive_model(rows: list[dict[str, object]]) -> dict[str, object]:
    y = np.asarray([float(row["log_resource_expansion"]) for row in rows], dtype=float)
    x = np.column_stack(
        [
            np.ones(len(rows)),
            np.log1p(
                np.asarray([float(row["host_family_count"]) for row in rows], dtype=float)
            ),
            np.log1p(
                np.asarray([int(row["native_resource_units"]) for row in rows], dtype=float)
            ),
            np.log1p(
                np.asarray([int(row["resolved_host_species"]) for row in rows], dtype=float)
            ),
        ]
    )
    coef, *_ = np.linalg.lstsq(x, y, rcond=None)
    fitted = x @ coef
    ss_total = float(np.sum((y - np.mean(y)) ** 2))
    ss_resid = float(np.sum((y - fitted) ** 2))
    r2 = None if ss_total <= 0 else 1.0 - ss_resid / ss_total
    return {
        "response": "log_resource_expansion",
        "predictors": [
            "intercept",
            "log1p(host_family_count)",
            "log1p(native_resource_units)",
            "log1p(resolved_host_species)",
        ],
        "coefficients": [float(value) for value in coef],
        "r_squared": None if r2 is None else float(r2),
        "inference": "descriptive least squares only",
    }


def summarize_rows(rows: list[dict[str, object]]) -> dict[str, object]:
    host = [float(row["host_family_count"]) for row in rows]
    added = [int(row["introduced_added_units"]) for row in rows]
    log_exp = [float(row["log_resource_expansion"]) for row in rows]
    resolved = [int(row["resolved_host_species"]) for row in rows]
    native = [int(row["native_resource_units"]) for row in rows]

    bands = {}
    for label in ("1_family", "2_families", "3_to_5_families", "6plus_families"):
        subset = [row for row in rows if row["host_breadth_stratum"] == label]
        if not subset:
            continue
        bands[label] = {
            "species": len(subset),
            "species_expanded": sum(
                int(row["introduced_added_units"] > 0) for row in subset
            ),
            "median_native_resource_units": float(
                np.median([row["native_resource_units"] for row in subset])
            ),
            "median_contemporary_resource_units": float(
                np.median([row["contemporary_resource_units"] for row in subset])
            ),
            "median_introduced_added_units": float(
                np.median([row["introduced_added_units"] for row in subset])
            ),
            "median_contemporary_over_native_ratio": float(
                np.median([row["contemporary_over_native_ratio"] for row in subset])
            ),
            "median_log_resource_expansion": float(
                np.median([row["log_resource_expansion"] for row in subset])
            ),
            "median_introduced_share_of_contemporary": float(
                np.median([row["introduced_share_of_contemporary"] for row in subset])
            ),
        }

    return {
        "species": len(rows),
        "species_expanded": sum(int(value > 0) for value in added),
        "fraction_species_expanded": (
            None if not rows else sum(int(value > 0) for value in added) / len(rows)
        ),
        "total_native_species_units": int(sum(native)),
        "total_contemporary_species_units": int(
            sum(int(row["contemporary_resource_units"]) for row in rows)
        ),
        "total_introduced_added_species_units": int(sum(added)),
        "spearman": {
            "host_family_count_vs_log_resource_expansion": spearman(host, log_exp),
            "host_family_count_vs_introduced_added_units": spearman(host, added),
            "resolved_host_species_vs_log_resource_expansion": spearman(
                resolved, log_exp
            ),
            "native_resource_units_vs_log_resource_expansion": spearman(
                native, log_exp
            ),
        },
        "host_breadth_strata": bands,
        "descriptive_model": descriptive_model(rows),
    }


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--protocol-json", type=Path, required=True)
    ap.add_argument("--descriptors-csv", type=Path, required=True)
    ap.add_argument("--native-footprints-json", type=Path, required=True)
    ap.add_argument("--insect-host-csv", type=Path, required=True)
    ap.add_argument("--contemporary-distribution-csv", type=Path, required=True)
    ap.add_argument("--output-csv", type=Path, required=True)
    ap.add_argument("--output-json", type=Path, required=True)
    args = ap.parse_args()

    protocol = json.loads(args.protocol_json.read_text(encoding="utf-8"))
    if protocol.get("schema") != (
        "ttf_butterfly_anthropogenic_resource_expansion_protocol_v0.1"
    ):
        raise RuntimeError("unexpected resource-expansion protocol schema")
    if protocol.get("status") != (
        "FROZEN_EXPLORATORY_BEFORE_FULL_S1_CONTEMPORARY_RESOURCE_RESULT"
    ):
        raise RuntimeError("resource-expansion protocol is not frozen")

    descriptors = load_descriptors(args.descriptors_csv)
    native = load_native_footprints(args.native_footprints_json)
    contemporary, diagnostics = load_contemporary_sidecars(
        args.insect_host_csv,
        args.contemporary_distribution_csv,
    )

    rows = []
    for name, descriptor in descriptors.items():
        host_families = float(descriptor["host_family_count"])
        native_units = native.get(name, frozenset())
        if host_families <= 0 or len(native_units) <= 0:
            continue
        contemporary_units = contemporary.get(name, frozenset())
        if not native_units.issubset(contemporary_units):
            raise RuntimeError(f"native footprint not subset of contemporary: {name}")

        n_native = len(native_units)
        n_contemporary = len(contemporary_units)
        added = n_contemporary - n_native
        resolved_hosts = int(descriptor["resolved_host_species"])
        rows.append(
            {
                "species": name,
                "host_family_count": host_families,
                "host_breadth_stratum": host_band(host_families),
                "resolved_host_species": resolved_hosts,
                "hosts_with_primary_native_units": int(
                    descriptor["hosts_with_primary_native_units"]
                ),
                "hosts_with_contemporary_units": int(
                    diagnostics.get(name, {}).get(
                        "hosts_with_primary_native_units", 0
                    )
                ),
                "native_resource_units": n_native,
                "contemporary_resource_units": n_contemporary,
                "introduced_added_units": added,
                "contemporary_over_native_ratio": n_contemporary / n_native,
                "log_resource_expansion": (
                    math.log1p(n_contemporary) - math.log1p(n_native)
                ),
                "introduced_share_of_contemporary": (
                    0.0 if n_contemporary == 0 else added / n_contemporary
                ),
                "host_taxonomy_lower_bound_adequate": (
                    resolved_hosts >= host_families
                ),
                "wing_size_proxy": float(descriptor["wing_size_proxy"]),
                "voltinism": str(descriptor["voltinism"]),
            }
        )

    rows = sorted(rows, key=lambda row: str(row["species"]))
    if len(rows) != 239:
        raise RuntimeError(
            f"expected 239 resource-eligible S1 species; reconstructed {len(rows)}"
        )

    adequate = [
        row
        for row in rows
        if bool(row["host_taxonomy_lower_bound_adequate"])
    ]

    payload = {
        "schema": "ttf_butterfly_anthropogenic_resource_expansion_result_v0.1",
        "status": "EXPLORATORY_FULL_S1_RESOURCE_RESULT",
        "protocol": str(args.protocol_json),
        "source_identity": {
            "descriptor_sha256": EXPECTED_DESCRIPTOR_SHA256,
            "native_footprint_sha256": EXPECTED_NATIVE_FOOTPRINT_SHA256,
        },
        "full_resource_eligible_panel": summarize_rows(rows),
        "host_taxonomy_lower_bound_adequate_subset": summarize_rows(adequate),
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
    print(json.dumps(payload, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
