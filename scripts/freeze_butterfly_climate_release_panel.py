#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from collections import Counter
from dataclasses import asdict
from pathlib import Path

from ttf.butterfly_climate_release import (
    ExpansionDescriptor,
    host_breadth_stratum,
    select_independent_panel,
    species_digest,
)


EXPECTED_DESCRIPTOR_SHA256 = (
    "894f48dbca1760fc4fa75bfee8f663540ab4b9380f8b9daf2bc09440e8bb0cdc"
)
EXPECTED_DESCRIPTOR_RUN_ID = 36222306369


def sha256_path(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def load_descriptors(path: Path) -> list[ExpansionDescriptor]:
    if sha256_path(path) != EXPECTED_DESCRIPTOR_SHA256:
        raise RuntimeError("S1 resource descriptor SHA drift")
    rows = []
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        required = {
            "species",
            "host_family_count",
            "host_wgsrpd3_unit_count",
            "resolved_host_species",
            "wing_size_proxy",
            "voltinism",
        }
        if not required <= set(reader.fieldnames or ()):
            raise RuntimeError("descriptor schema drift")
        for row in reader:
            rows.append(
                ExpansionDescriptor(
                    species=str(row["species"]).strip(),
                    host_family_count=float(row["host_family_count"]),
                    host_wgsrpd3_unit_count=int(row["host_wgsrpd3_unit_count"]),
                    resolved_host_species=int(row["resolved_host_species"]),
                    wing_size_proxy=float(row["wing_size_proxy"]),
                    voltinism=str(row["voltinism"]).strip(),
                )
            )
    if len(rows) != 339 or len({row.species for row in rows}) != 339:
        raise RuntimeError("expected exact 339-species S1 descriptor table")
    return rows


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--descriptors-csv", type=Path, required=True)
    ap.add_argument("--protocol-json", type=Path, required=True)
    ap.add_argument("--output-json", type=Path, required=True)
    ap.add_argument("--output-csv", type=Path, required=True)
    args = ap.parse_args()

    protocol = json.loads(args.protocol_json.read_text(encoding="utf-8"))
    if protocol.get("schema") != "ttf_butterfly_climate_release_independent_test_v0.1":
        raise RuntimeError("unexpected independent-test protocol schema")
    if protocol.get("status") != (
        "FROZEN_PILOT_DERIVED_HYPOTHESIS_BEFORE_INDEPENDENT_GBIF_OR_CLIMATE"
    ):
        raise RuntimeError("independent-test protocol is not frozen")

    descriptors = load_descriptors(args.descriptors_csv)
    excluded = tuple(map(str, protocol["independent_panel"]["excluded_species"]))
    target_total = int(protocol["independent_panel"]["target_species_total"])
    strata_rule = protocol["independent_panel"]["host_breadth_strata"]
    per_stratum_values = {int(v["target_species"]) for v in strata_rule.values()}
    if len(per_stratum_values) != 1:
        raise RuntimeError("independent panel requires equal stratum targets")
    per_stratum = next(iter(per_stratum_values))
    minimum_native_units = 10

    panel = select_independent_panel(
        descriptors,
        excluded_species=excluded,
        per_stratum=per_stratum,
        minimum_native_units=minimum_native_units,
    )
    if len(panel) != target_total:
        raise RuntimeError("selected panel size drift")
    if set(panel) & set(excluded):
        raise RuntimeError("pilot species leaked into independent panel")

    by_species = {row.species: row for row in descriptors}
    selected = [by_species[name] for name in panel]
    stratum_counts = Counter(
        host_breadth_stratum(row.host_family_count)
        for row in selected
    )
    expected_counts = {
        name: int(rule["target_species"])
        for name, rule in strata_rule.items()
    }
    if dict(sorted(stratum_counts.items())) != dict(sorted(expected_counts.items())):
        raise RuntimeError("host-breadth stratum balance drift")

    payload = {
        "schema": "ttf_butterfly_climate_release_independent_panel_v0.1",
        "status": "FROZEN_BEFORE_INDEPENDENT_GBIF_OR_CLIMATE",
        "protocol": str(args.protocol_json),
        "source": {
            "descriptor_workflow_run_id": EXPECTED_DESCRIPTOR_RUN_ID,
            "descriptor_sha256": EXPECTED_DESCRIPTOR_SHA256,
            "descriptor_species": len(descriptors),
        },
        "selection": {
            "pilot_species_excluded": list(excluded),
            "minimum_native_host_resource_units": minimum_native_units,
            "host_taxonomy_lower_bound_required": True,
            "per_stratum": per_stratum,
            "stratum_counts": dict(sorted(stratum_counts.items())),
            "new_butterfly_gbif_used": False,
            "climate_response_used": False,
            "genetic_response_used": False,
        },
        "species_count": len(panel),
        "species": list(panel),
        "species_sha256": species_digest(panel),
        "descriptors": [
            {
                **asdict(row),
                "host_breadth_stratum": host_breadth_stratum(
                    row.host_family_count
                ),
            }
            for row in selected
        ],
    }
    args.output_json.parent.mkdir(parents=True, exist_ok=True)
    args.output_json.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )

    fields = [
        "panel_order",
        "species",
        "host_breadth_stratum",
        "host_family_count",
        "host_wgsrpd3_unit_count",
        "resolved_host_species",
        "wing_size_proxy",
        "voltinism",
    ]
    with args.output_csv.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for i, row in enumerate(selected):
            writer.writerow(
                {
                    "panel_order": i,
                    "species": row.species,
                    "host_breadth_stratum": host_breadth_stratum(
                        row.host_family_count
                    ),
                    "host_family_count": row.host_family_count,
                    "host_wgsrpd3_unit_count": row.host_wgsrpd3_unit_count,
                    "resolved_host_species": row.resolved_host_species,
                    "wing_size_proxy": row.wing_size_proxy,
                    "voltinism": row.voltinism,
                }
            )

    print(
        json.dumps(
            {
                "species_count": len(panel),
                "species_sha256": payload["species_sha256"],
                "stratum_counts": payload["selection"]["stratum_counts"],
                "species": list(panel),
            },
            indent=2,
            sort_keys=True,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
