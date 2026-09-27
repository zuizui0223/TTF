#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from collections import defaultdict
from pathlib import Path

import numpy as np


EXPECTED_MECHANISM_RESULT_SHA256 = (
    "37fb75f7ff500528a6e68a05ac8088acf4c1b77687577a757a77e97da05bc3a2"
)
EXPECTED_MECHANISM_METRICS_SHA256 = (
    "b0f16c5fa9a5b4a0842d6d23f69de7a1f5e938a4a96fea426c97df2dd73e63aa"
)


def sha256_path(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def require_sha(path: Path, expected: str, label: str) -> None:
    observed = sha256_path(path)
    if observed != expected:
        raise RuntimeError(
            f"{label} SHA drift: expected {expected}, observed {observed}"
        )


def load_host_metadata(path: Path) -> dict[str, dict[str, str]]:
    result: dict[str, dict[str, str]] = {}
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        required = {"accepted_plant_name_id", "accepted_name", "family"}
        if not required <= set(reader.fieldnames or ()):
            raise RuntimeError("insect-host metadata sidecar schema drift")
        for row in reader:
            host_id = str(row["accepted_plant_name_id"]).strip()
            if not host_id:
                continue
            value = {
                "accepted_name": str(row.get("accepted_name", "")).strip(),
                "family": str(row.get("family", "")).strip(),
            }
            if host_id in result and result[host_id] != value:
                raise RuntimeError(f"host metadata conflict for {host_id}")
            result[host_id] = value
    return result


def load_units(path: Path) -> dict[str, set[str]]:
    result: dict[str, set[str]] = defaultdict(set)
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        required = {"accepted_plant_name_id", "area_code_l3"}
        if not required <= set(reader.fieldnames or ()):
            raise RuntimeError("host distribution sidecar schema drift")
        for row in reader:
            host_id = str(row["accepted_plant_name_id"]).strip()
            unit = str(row["area_code_l3"]).strip()
            if host_id and unit:
                result[host_id].add(unit)
    return result


def load_species_metrics(path: Path) -> dict[str, dict[str, object]]:
    require_sha(
        path,
        EXPECTED_MECHANISM_METRICS_SHA256,
        "mechanism species metrics",
    )
    result = {}
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        required = {
            "species",
            "host_family_count",
            "resolved_host_species",
            "introduced_added_units",
            "host_taxonomy_lower_bound_adequate",
        }
        if not required <= set(reader.fieldnames or ()):
            raise RuntimeError("mechanism metrics schema drift")
        for row in reader:
            name = str(row["species"]).strip()
            result[name] = {
                "species": name,
                "host_family_count": float(row["host_family_count"]),
                "resolved_host_species": int(row["resolved_host_species"]),
                "introduced_added_units": int(row["introduced_added_units"]),
                "host_taxonomy_lower_bound_adequate": (
                    str(row["host_taxonomy_lower_bound_adequate"]).strip()
                    == "True"
                ),
            }
    if len(result) != 239:
        raise RuntimeError(f"expected 239 mechanism species; got {len(result)}")
    return result


def concentration_summary(
    plant_rows: list[dict[str, object]],
) -> dict[str, object]:
    credits = np.asarray(
        [float(row["total_fractional_butterfly_unit_credit"]) for row in plant_rows],
        dtype=float,
    )
    total = float(np.sum(credits))
    if total <= 0:
        raise RuntimeError("plant credit total is non-positive")
    ordered = np.sort(credits)[::-1]
    cumulative = np.cumsum(ordered) / total

    def share_top(n: int) -> float:
        return float(np.sum(ordered[:n]) / total)

    def hosts_to_reach(target: float) -> int:
        return int(np.searchsorted(cumulative, target, side="left") + 1)

    return {
        "total_fractional_butterfly_unit_credit": total,
        "host_plants_with_positive_credit": int(len(plant_rows)),
        "share_top_1_host": share_top(1),
        "share_top_5_hosts": share_top(5),
        "share_top_10_hosts": share_top(10),
        "share_top_25_hosts": share_top(25),
        "hosts_to_reach_25_percent": hosts_to_reach(0.25),
        "hosts_to_reach_50_percent": hosts_to_reach(0.50),
        "hosts_to_reach_75_percent": hosts_to_reach(0.75),
    }


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--protocol-json", type=Path, required=True)
    ap.add_argument("--mechanism-result-json", type=Path, required=True)
    ap.add_argument("--species-metrics-csv", type=Path, required=True)
    ap.add_argument("--insect-host-csv", type=Path, required=True)
    ap.add_argument("--native-distribution-csv", type=Path, required=True)
    ap.add_argument("--contemporary-distribution-csv", type=Path, required=True)
    ap.add_argument("--output-plant-csv", type=Path, required=True)
    ap.add_argument("--output-butterfly-examples-csv", type=Path, required=True)
    ap.add_argument("--output-json", type=Path, required=True)
    args = ap.parse_args()

    protocol = json.loads(args.protocol_json.read_text(encoding="utf-8"))
    if protocol.get("schema") != "ttf_butterfly_host_plant_driver_analysis_v0.1":
        raise RuntimeError("unexpected plant-driver protocol schema")
    if protocol.get("status") != (
        "FROZEN_POST_RESULT_EXPLORATORY_BEFORE_PLANT_LEVEL_AGGREGATION"
    ):
        raise RuntimeError("plant-driver protocol is not frozen")

    require_sha(
        args.mechanism_result_json,
        EXPECTED_MECHANISM_RESULT_SHA256,
        "mechanism result",
    )
    mechanism = json.loads(args.mechanism_result_json.read_text(encoding="utf-8"))
    if mechanism.get("schema") != (
        "ttf_butterfly_resource_expansion_mechanism_result_v0.1"
    ):
        raise RuntimeError("unexpected mechanism result schema")

    credits_by_butterfly = {
        str(species): {
            str(host_id): float(value)
            for host_id, value in credits.items()
        }
        for species, credits in mechanism.get(
            "host_fractional_credits_for_expanded_species", {}
        ).items()
    }
    if len(credits_by_butterfly) != 206:
        raise RuntimeError(
            f"expected 206 expanded butterflies with credits; "
            f"got {len(credits_by_butterfly)}"
        )

    species_metrics = load_species_metrics(args.species_metrics_csv)
    host_meta = load_host_metadata(args.insect_host_csv)
    native = load_units(args.native_distribution_csv)
    contemporary = load_units(args.contemporary_distribution_csv)

    aggregate: dict[str, dict[str, object]] = {}
    dominant_host_by_butterfly: dict[str, str] = {}

    for butterfly, credits in sorted(credits_by_butterfly.items()):
        added_units = int(species_metrics[butterfly]["introduced_added_units"])
        if added_units <= 0:
            raise RuntimeError(f"positive-credit butterfly has no added units: {butterfly}")
        if abs(sum(credits.values()) - added_units) > 1e-8:
            raise RuntimeError(f"credit sum drift for {butterfly}")

        maximum = max(credits.values())
        top_hosts = sorted(
            host_id
            for host_id, credit in credits.items()
            if abs(credit - maximum) <= 1e-12
        )
        dominant_host_by_butterfly[butterfly] = top_hosts[0]

        for host_id, credit in credits.items():
            if host_id not in aggregate:
                meta = host_meta.get(
                    host_id,
                    {"accepted_name": "", "family": ""},
                )
                native_units = set(native.get(host_id, set()))
                contemporary_units = set(contemporary.get(host_id, set()))
                if not native_units.issubset(contemporary_units):
                    raise RuntimeError(
                        f"native not subset contemporary for host {host_id}"
                    )
                aggregate[host_id] = {
                    "accepted_plant_name_id": host_id,
                    "accepted_name": meta["accepted_name"],
                    "family": meta["family"],
                    "butterflies": set(),
                    "shares": [],
                    "total_fractional_butterfly_unit_credit": 0.0,
                    "butterfly_species_where_host_is_top_contributor_including_ties": 0,
                    "butterfly_species_where_host_is_sole_contributor": 0,
                    "native_wgsrpd3_units": len(native_units),
                    "contemporary_wgsrpd3_units": len(contemporary_units),
                    "introduced_added_wgsrpd3_units": len(
                        contemporary_units - native_units
                    ),
                }

            row = aggregate[host_id]
            row["butterflies"].add(butterfly)
            row["shares"].append(float(credit) / added_units)
            row["total_fractional_butterfly_unit_credit"] += float(credit)
            if host_id in top_hosts:
                row[
                    "butterfly_species_where_host_is_top_contributor_including_ties"
                ] += 1
            if len(credits) == 1:
                row["butterfly_species_where_host_is_sole_contributor"] += 1

    plant_rows = []
    for host_id, row in aggregate.items():
        shares = list(row.pop("shares"))
        butterflies = sorted(row.pop("butterflies"))
        plant_rows.append(
            {
                **row,
                "butterfly_species_with_positive_fractional_credit": len(butterflies),
                "median_fractional_share_within_butterfly": float(
                    np.median(np.asarray(shares, dtype=float))
                ),
                "maximum_fractional_share_within_butterfly": float(max(shares)),
                "butterfly_examples": "; ".join(butterflies[:12]),
            }
        )

    plant_rows = sorted(
        plant_rows,
        key=lambda row: (
            -int(row["butterfly_species_with_positive_fractional_credit"]),
            -float(row["total_fractional_butterfly_unit_credit"]),
            str(row["accepted_name"]),
            str(row["accepted_plant_name_id"]),
        ),
    )

    overall_by_credit = sorted(
        plant_rows,
        key=lambda row: (
            -float(row["total_fractional_butterfly_unit_credit"]),
            -int(row["butterfly_species_with_positive_fractional_credit"]),
            str(row["accepted_name"]),
        ),
    )

    def example_row(butterfly: str) -> dict[str, object]:
        metrics = species_metrics[butterfly]
        dominant_id = dominant_host_by_butterfly[butterfly]
        credits = credits_by_butterfly[butterfly]
        added_units = int(metrics["introduced_added_units"])
        meta = host_meta.get(dominant_id, {"accepted_name": "", "family": ""})
        return {
            "butterfly_species": butterfly,
            "host_family_count": metrics["host_family_count"],
            "resolved_host_species": metrics["resolved_host_species"],
            "introduced_added_resource_units": added_units,
            "dominant_host_accepted_name": meta["accepted_name"],
            "dominant_host_family": meta["family"],
            "dominant_host_fractional_share": (
                float(credits[dominant_id]) / added_units
            ),
            "dominant_host_introduced_added_wgsrpd3_units": len(
                set(contemporary.get(dominant_id, set()))
                - set(native.get(dominant_id, set()))
            ),
        }

    adequate_expanded = [
        metrics
        for metrics in species_metrics.values()
        if metrics["host_taxonomy_lower_bound_adequate"]
        and int(metrics["introduced_added_units"]) > 0
    ]
    one_family = sorted(
        (
            metrics
            for metrics in adequate_expanded
            if float(metrics["host_family_count"]) == 1.0
        ),
        key=lambda row: (
            -int(row["introduced_added_units"]),
            str(row["species"]),
        ),
    )[:10]
    top_all = sorted(
        adequate_expanded,
        key=lambda row: (
            -int(row["introduced_added_units"]),
            str(row["species"]),
        ),
    )[:10]

    butterfly_examples = []
    for group, rows in (
        ("one_family_top_expansion", one_family),
        ("overall_top_expansion", top_all),
    ):
        for rank, metrics in enumerate(rows, start=1):
            butterfly_examples.append(
                {
                    "example_group": group,
                    "rank": rank,
                    **example_row(str(metrics["species"])),
                }
            )

    concentration = concentration_summary(plant_rows)
    top_n = int(protocol["rankings"]["top_rows_to_report"])
    payload = {
        "schema": "ttf_butterfly_host_plant_driver_result_v0.1",
        "status": "EXPLORATORY_PLANT_LEVEL_DRIVER_AGGREGATION_COMPLETE",
        "protocol": str(args.protocol_json),
        "species": {
            "resource_eligible_butterflies": 239,
            "expanded_butterflies": 206,
        },
        "host_plants": {
            "positive_credit_host_plants": len(plant_rows),
            "concentration": concentration,
            "top_by_butterfly_reach": plant_rows[:top_n],
            "top_by_total_fractional_credit": overall_by_credit[:top_n],
        },
        "butterfly_examples": butterfly_examples,
        "ecological_interpretation": {
            "question": (
                "Whether reconstructed anthropogenic butterfly resource "
                "opportunity is broadly distributed among host plants or "
                "concentrated in a small set of widely redistributed plants."
            ),
            "claim_boundary": protocol["interpretation_boundary"],
        },
    }

    args.output_plant_csv.parent.mkdir(parents=True, exist_ok=True)
    plant_fields = list(plant_rows[0].keys())
    with args.output_plant_csv.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=plant_fields)
        writer.writeheader()
        writer.writerows(plant_rows)

    example_fields = list(butterfly_examples[0].keys())
    with args.output_butterfly_examples_csv.open(
        "w", newline="", encoding="utf-8"
    ) as handle:
        writer = csv.DictWriter(handle, fieldnames=example_fields)
        writer.writeheader()
        writer.writerows(butterfly_examples)

    args.output_json.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    print(
        json.dumps(
            {
                "schema": payload["schema"],
                "status": payload["status"],
                "positive_credit_host_plants": len(plant_rows),
                "concentration": concentration,
                "top_by_butterfly_reach": [
                    {
                        key: row[key]
                        for key in (
                            "accepted_name",
                            "family",
                            "butterfly_species_with_positive_fractional_credit",
                            "total_fractional_butterfly_unit_credit",
                            "introduced_added_wgsrpd3_units",
                        )
                    }
                    for row in plant_rows[:10]
                ],
            },
            indent=2,
            sort_keys=True,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
