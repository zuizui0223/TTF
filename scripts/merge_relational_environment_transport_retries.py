#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import json
from collections import Counter
from pathlib import Path


FIELDS = ["species", "source_key", "latitude", "longitude", "priority_rank", "priority_sha256"]


def load_base(input_dir: Path):
    species = {}
    rows_by_species = {}
    for path in sorted(input_dir.glob("ledger-*.json")):
        payload = json.loads(path.read_text())
        if payload.get("schema") != "ttf_relational_environment_occurrence_shard_v0.2":
            raise RuntimeError(f"unexpected base shard schema: {path}")
        for row in payload["species"]:
            name = str(row["species"])
            if name in species:
                raise RuntimeError(f"duplicate base species: {name}")
            species[name] = row
    for path in sorted(input_dir.glob("occurrences-*.csv")):
        for row in csv.DictReader(path.open(encoding="utf-8")):
            rows_by_species.setdefault(str(row["species"]), []).append(dict(row))
    if len(species) != 1000:
        raise RuntimeError(f"expected 1000 base species, found {len(species)}")
    return species, rows_by_species


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--base-dir", type=Path, required=True)
    ap.add_argument("--retry-ledger", type=Path, action="append", default=[])
    ap.add_argument("--retry-csv", type=Path, action="append", default=[])
    ap.add_argument("--output-csv", type=Path, required=True)
    ap.add_argument("--output-ledger", type=Path, required=True)
    args = ap.parse_args()

    species, rows_by_species = load_base(args.base_dir)
    for ledger_path, csv_path in zip(args.retry_ledger, args.retry_csv):
        payload = json.loads(ledger_path.read_text())
        if payload.get("schema") != "ttf_relational_environment_transport_retry_v0.1":
            raise RuntimeError("unexpected retry schema")
        if any(bool(v) for v in payload["response_firewall"].values()):
            raise RuntimeError("retry firewall open")
        retry_rows = {}
        for row in csv.DictReader(csv_path.open(encoding="utf-8")):
            retry_rows.setdefault(str(row["species"]), []).append(dict(row))
        for row in payload["species"]:
            name = str(row["species"])
            if species[name]["status"] != "REQUEST_ERROR":
                continue
            species[name] = row
            rows_by_species[name] = retry_rows.get(name, [])

    unresolved = [name for name, row in species.items() if row["status"] == "REQUEST_ERROR"]
    all_rows = []
    for name in sorted(rows_by_species):
        all_rows.extend(rows_by_species[name])
    all_rows.sort(key=lambda row: (row["species"], int(row["priority_rank"]), int(row["source_key"])))

    args.output_csv.parent.mkdir(parents=True, exist_ok=True)
    with args.output_csv.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDS)
        writer.writeheader()
        writer.writerows(all_rows)

    counts = Counter(str(row["status"]) for row in species.values())
    passed = counts.get("PASS_OCCURRENCE_GEOMETRY", 0)
    payload = {
        "schema": "ttf_relational_environment_occurrence_acquisition_v0.2",
        "status": (
            "INCOMPLETE_TECHNICAL_EXECUTION"
            if unresolved else
            ("PASS_TO_CHELSA_EXTRACTION" if passed >= 500 else "NOT_EVALUABLE_OCCURRENCE_GEOMETRY")
        ),
        "species": 1000,
        "species_passing_ge_30": passed,
        "status_counts": dict(sorted(counts.items())),
        "retained_rows": len(all_rows),
        "unresolved_request_error_species": unresolved,
        "response_firewall": {
            "study_B_sequence_identity_opened": False,
            "study_B_pairwise_genetic_distances_opened": False,
            "study_B_T_st_computed": False,
            "study_B_beta_R_computed": False,
        },
    }
    args.output_ledger.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    print(json.dumps({
        "status": payload["status"],
        "species_passing_ge_30": passed,
        "request_errors": len(unresolved),
    }, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
