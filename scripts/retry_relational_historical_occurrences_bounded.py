#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

try:
    from scripts.acquire_relational_historical_occurrences_bounded import FIELDS, bounded_fetch
except ModuleNotFoundError:
    from acquire_relational_historical_occurrences_bounded import FIELDS, bounded_fetch


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--plan", type=Path, required=True)
    ap.add_argument("--output-csv", type=Path, required=True)
    ap.add_argument("--output-ledger", type=Path, required=True)
    ap.add_argument("--batch-index", type=int, required=True)
    ap.add_argument("--batch-size", type=int, default=4)
    ap.add_argument("--per-species-timeout-seconds", type=int, default=2700)
    args = ap.parse_args()

    plan = json.loads(args.plan.read_text())
    if plan.get("schema") != "ttf_relational_historical_occurrence_retry_plan_v0.2":
        raise RuntimeError("unexpected Study-C bounded retry plan")
    if plan.get("status") not in {
        "RETRY_EXACT_INITIAL_REQUEST_ERRORS",
        "NO_RETRY_NEEDED",
    }:
        raise RuntimeError("unexpected Study-C bounded retry plan status")
    if plan.get("scientific_query_change") is not False:
        raise RuntimeError("Study-C retry plan changed scientific query")
    if any(bool(v) for v in plan["response_firewall"].values()):
        raise RuntimeError("Study-C retry-plan response firewall is open")

    names = list(map(str, plan["request_error_species"]))
    if args.batch_size != 4:
        raise RuntimeError("Study-C bounded retry requires exact batch size 4")
    if args.per_species_timeout_seconds != 2700:
        raise RuntimeError("Study-C bounded retry requires exact 2700-second per-species timeout")

    start = args.batch_index * args.batch_size
    batch = names[start : start + args.batch_size]
    if not batch:
        raise RuntimeError(f"empty retry batch {args.batch_index}")

    retained_rows: list[dict[str, object]] = []
    ledgers: list[dict[str, object]] = []
    for ordinal, species in enumerate(batch):
        retained, ledger = bounded_fetch(species, args.per_species_timeout_seconds)
        retained_rows.extend(retained)
        ledger = dict(ledger)
        ledger["bounded_retry_batch_index"] = int(args.batch_index)
        ledger["bounded_retry_batch_ordinal"] = int(ordinal)
        ledger["per_species_wall_timeout_seconds"] = int(args.per_species_timeout_seconds)
        ledgers.append(ledger)
        print(json.dumps({
            "retry_batch_index": args.batch_index,
            "ordinal": ordinal,
            "species": species,
            "status": ledger["status"],
        }, sort_keys=True), flush=True)

    args.output_csv.parent.mkdir(parents=True, exist_ok=True)
    with args.output_csv.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDS)
        writer.writeheader()
        writer.writerows(retained_rows)

    payload = {
        "schema": "ttf_relational_historical_occurrence_bounded_retry_batch_v0.2",
        "retry_round": 1,
        "batch_index": int(args.batch_index),
        "batch_size": int(args.batch_size),
        "species_requested": batch,
        "species": ledgers,
        "retained_occurrence_rows": len(retained_rows),
        "scientific_query_change": False,
        "transport_only_change": {
            "per_species_wall_timeout_seconds": int(args.per_species_timeout_seconds),
            "salvage_completed_species_even_when_sibling_times_out": True,
        },
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
