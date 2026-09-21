#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

try:
    from scripts.acquire_relational_environment_occurrences import fetch_species
except ModuleNotFoundError:
    from acquire_relational_environment_occurrences import fetch_species

from ttf.relational_external import shard_items


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--candidates", type=Path, required=True)
    ap.add_argument("--output-csv", type=Path, required=True)
    ap.add_argument("--output-ledger", type=Path, required=True)
    ap.add_argument("--shard-index", type=int, default=0)
    ap.add_argument("--shards", type=int, default=1)
    args = ap.parse_args()

    rows = list(csv.DictReader(args.candidates.open(encoding="utf-8")))
    names = [str(row["species"]).strip() for row in rows]
    if not (500 <= len(names) <= 1000) or len(set(names)) != len(names):
        raise RuntimeError("Study C candidates must be 500-1000 unique frozen species")
    shard = list(shard_items(names, shard_index=args.shard_index, shards=args.shards))

    retained_rows = []
    ledgers = []
    for i, species in enumerate(shard, start=1):
        try:
            retained, ledger = fetch_species(species)
        except Exception as exc:
            retained, ledger = [], {
                "species": species,
                "status": "REQUEST_ERROR",
                "error": f"{type(exc).__name__}: {exc}",
            }
        retained_rows.extend(retained)
        ledgers.append(ledger)
        if i % 5 == 0:
            print(json.dumps({"shard": args.shard_index, "completed": i, "total": len(shard)}))

    args.output_csv.parent.mkdir(parents=True, exist_ok=True)
    fields = ["species", "source_key", "latitude", "longitude", "priority_rank", "priority_sha256"]
    with args.output_csv.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(retained_rows)

    payload = {
        "schema": "ttf_relational_historical_occurrence_shard_v0.1",
        "shard_index": args.shard_index,
        "shards": args.shards,
        "candidate_species_total": len(names),
        "species_requested": len(shard),
        "retained_occurrence_rows": len(retained_rows),
        "species": ledgers,
        "inheritance": "Exact Study-B GBIF acquisition function reused without ecological/genetic-response changes.",
        "response_firewall": {
            "Study_C_sequence_identity_opened": False,
            "Study_C_pairwise_genetic_distances_opened": False,
            "Study_C_T_st_computed": False,
            "Study_C_beta_hist_computed": False,
        },
    }
    args.output_ledger.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
