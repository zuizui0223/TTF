#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import json
from collections import defaultdict
from pathlib import Path


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--native-unit-table", type=Path, required=True)
    ap.add_argument("--contemporary-footprints-json", type=Path, required=True)
    ap.add_argument("--output-json", type=Path, required=True)
    args = ap.parse_args()

    contemporary_payload = json.loads(
        args.contemporary_footprints_json.read_text(encoding="utf-8")
    )
    if contemporary_payload.get("schema") != (
        "ttf_butterfly_resource_envelope_contemporary_footprints_v0.1"
    ):
        raise RuntimeError("unexpected contemporary footprint schema")
    contemporary = {
        str(name): set(map(str, units))
        for name, units in contemporary_payload.get("species", {}).items()
    }

    native = defaultdict(set)
    observed = defaultdict(set)
    record_count = defaultdict(int)
    with args.native_unit_table.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        required = {
            "species",
            "wgsrpd3_code",
            "host_available",
            "butterfly_observed",
            "butterfly_record_count",
        }
        if not required <= set(reader.fieldnames or ()):
            raise RuntimeError("native unit table schema drift")
        for row in reader:
            species = str(row["species"]).strip()
            code = str(row["wgsrpd3_code"]).strip()
            if int(row["host_available"]):
                native[species].add(code)
            if int(row["butterfly_observed"]):
                observed[species].add(code)
            record_count[species] += int(row["butterfly_record_count"])

    species = sorted(set(native) | set(observed) | set(contemporary))
    rows = []
    for name in species:
        n = native[name]
        c = contemporary.get(name, set())
        o = observed[name]
        if not n.issubset(c):
            raise RuntimeError(f"native footprint not subset of contemporary: {name}")
        within_native = o & n
        outside_native = o - n
        within_contemporary = o & c
        outside_contemporary = o - c
        explained_by_introduced = outside_native & c
        rows.append(
            {
                "species": name,
                "butterfly_observed_units": len(o),
                "native_host_units": len(n),
                "contemporary_host_units": len(c),
                "observed_within_native": len(within_native),
                "observed_outside_native": len(outside_native),
                "observed_within_contemporary": len(within_contemporary),
                "observed_outside_contemporary": len(outside_contemporary),
                "outside_native_units_explained_by_introduced_host_ranges": len(
                    explained_by_introduced
                ),
                "fraction_outside_native_explained_by_introduced": (
                    None
                    if not outside_native
                    else len(explained_by_introduced) / len(outside_native)
                ),
                "occurrence_records": int(record_count[name]),
            }
        )

    total_outside_native = sum(row["observed_outside_native"] for row in rows)
    total_explained = sum(
        row["outside_native_units_explained_by_introduced_host_ranges"]
        for row in rows
    )
    total_outside_contemporary = sum(
        row["observed_outside_contemporary"] for row in rows
    )
    payload = {
        "schema": "ttf_butterfly_resource_envelope_native_contemporary_occurrence_overlap_v0.1",
        "status": "EXPLORATORY_RESOURCE_ENVELOPE_VALIDITY_DIAGNOSTIC",
        "species": rows,
        "summary": {
            "species": len(rows),
            "observed_outside_native_species_units": total_outside_native,
            "outside_native_explained_by_introduced_host_ranges": total_explained,
            "fraction_outside_native_explained_by_introduced": (
                None
                if total_outside_native == 0
                else total_explained / total_outside_native
            ),
            "observed_outside_contemporary_species_units": total_outside_contemporary,
        },
        "interpretation": (
            "Observed butterfly units recovered only after introduced host ranges "
            "are admitted demonstrate that native-only host geography is not an "
            "appropriate contemporary resource envelope for modern GBIF records. "
            "Residual outside-contemporary units remain a diagnostic for host-record "
            "incompleteness, taxonomic mismatch, dispersal/vagrancy, mapping error, "
            "or resource use not represented in the frozen HOSTS-WCVP reconstruction."
        ),
        "claim_boundary": {
            "exploratory_only": True,
            "outside_contemporary_is_not_evidence_of_host_independence": True,
            "gbif_absence_is_not_true_absence": True,
            "graphium_may_be_missing_under_pre_resolver_repair_run": True,
        },
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
