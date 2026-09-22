#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from collections import Counter
from pathlib import Path


def sha256_path(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--input-dir", type=Path, required=True)
    ap.add_argument("--output-csv", type=Path, required=True)
    ap.add_argument("--output-ledger", type=Path, required=True)
    args = ap.parse_args()

    csv_paths = sorted(args.input_dir.glob("occurrences-*.csv"))
    ledger_paths = sorted(args.input_dir.glob("ledger-*.json"))
    if not csv_paths or not ledger_paths:
        raise RuntimeError("missing Study-C occurrence shards")

    rows = []
    for path in csv_paths:
        rows.extend(csv.DictReader(path.open(encoding="utf-8")))
    species_ledgers = []
    totals = set()
    for path in ledger_paths:
        payload = json.loads(path.read_text())
        if payload.get("schema") != "ttf_relational_historical_occurrence_shard_v0.1":
            raise RuntimeError(f"unexpected historical occurrence shard: {path}")
        if any(bool(v) for v in payload["response_firewall"].values()):
            raise RuntimeError("Study C response firewall is open")
        totals.add(int(payload["candidate_species_total"]))
        species_ledgers.extend(payload["species"])
    if len(totals) != 1:
        raise RuntimeError("candidate total drift across Study-C shards")
    candidate_total = next(iter(totals))
    names = [str(row["species"]) for row in species_ledgers]
    if len(names) != candidate_total or len(set(names)) != candidate_total:
        raise RuntimeError("Study-C occurrence acquisition did not cover exact candidate species")

    row_keys = [(row["species"], int(row["source_key"])) for row in rows]
    if len(row_keys) != len(set(row_keys)):
        raise RuntimeError("duplicate species/source_key rows across Study-C shards")
    rows.sort(key=lambda row: (row["species"], int(row["priority_rank"]), int(row["source_key"])))

    args.output_csv.parent.mkdir(parents=True, exist_ok=True)
    fields = ["species", "source_key", "latitude", "longitude", "priority_rank", "priority_sha256"]
    with args.output_csv.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)

    counts = Counter(str(row["status"]) for row in species_ledgers)
    passed = sum(row["status"] == "PASS_OCCURRENCE_GEOMETRY" for row in species_ledgers)
    payload = {
        "schema": "ttf_relational_historical_occurrence_acquisition_v0.1",
        "status": "PASS_TO_HISTORICAL_ASSET_EXTRACTION" if passed >= 500 else "NOT_EVALUABLE_HISTORICAL_OCCURRENCE_GEOMETRY",
        "candidate_species": candidate_total,
        "species_passing_ge_30": passed,
        "status_counts": dict(sorted(counts.items())),
        "retained_rows": len(rows),
        "occurrence_csv_sha256": sha256_path(args.output_csv),
        "response_firewall": {
            "Study_C_sequence_identity_opened": False,
            "Study_C_pairwise_genetic_distances_opened": False,
            "Study_C_T_st_computed": False,
            "Study_C_beta_hist_computed": False,
        },
    }
    args.output_ledger.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    print(json.dumps(payload, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
