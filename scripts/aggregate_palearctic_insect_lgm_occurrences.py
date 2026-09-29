#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--candidates", type=Path, required=True)
    ap.add_argument("--input-dir", type=Path, required=True)
    ap.add_argument("--output-csv", type=Path, required=True)
    ap.add_argument("--output-ledger", type=Path, required=True)
    args = ap.parse_args()

    candidate = json.loads(args.candidates.read_text())
    expected = list(map(str, candidate["species"]))
    expected_set = set(expected)

    ledger_rows: dict[str, dict] = {}
    for path in sorted(args.input_dir.rglob("ledger_*.json")):
        payload = json.loads(path.read_text())
        if payload.get("schema") != "ttf_palearctic_insect_lgm_gbif_shard_v0.1":
            continue
        for row in payload["species"]:
            species = str(row["species"])
            if species in ledger_rows:
                raise RuntimeError(f"duplicate ledger species: {species}")
            ledger_rows[species] = row

    if set(ledger_rows) != expected_set:
        missing = sorted(expected_set - set(ledger_rows))
        extra = sorted(set(ledger_rows) - expected_set)
        raise RuntimeError(f"GBIF shard coverage drift; missing={missing[:10]}, extra={extra[:10]}")

    occurrence_rows: list[dict[str, str]] = []
    for path in sorted(args.input_dir.rglob("occurrences_*.csv")):
        with path.open(newline="", encoding="utf-8") as handle:
            occurrence_rows.extend(csv.DictReader(handle))
    if any(row["species"] not in expected_set for row in occurrence_rows):
        raise RuntimeError("occurrence species outside frozen 143-species parent")

    occurrence_rows.sort(
        key=lambda row: (
            row["species"],
            int(row["priority_rank"]),
            int(row["source_key"]),
        )
    )
    args.output_csv.parent.mkdir(parents=True, exist_ok=True)
    fields = ["species", "source_key", "latitude", "longitude", "priority_rank", "priority_sha256"]
    with args.output_csv.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(occurrence_rows)

    status_counts: dict[str, int] = {}
    for row in ledger_rows.values():
        status = str(row["status"])
        status_counts[status] = status_counts.get(status, 0) + 1

    payload = {
        "schema": "ttf_palearctic_insect_lgm_gbif_aggregate_v0.1",
        "status": "COMPLETE_RESPONSE_BLIND_GBIF_CENSUS",
        "candidate_species": len(expected),
        "status_counts": dict(sorted(status_counts.items())),
        "retained_occurrence_rows": len(occurrence_rows),
        "species": [ledger_rows[name] for name in expected],
        "response_firewall": candidate["response_firewall"],
    }
    args.output_ledger.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
