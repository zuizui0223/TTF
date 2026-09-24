#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import json
import multiprocessing as mp
import queue
from pathlib import Path

try:
    from scripts.acquire_relational_environment_occurrences import fetch_species
except ModuleNotFoundError:
    from acquire_relational_environment_occurrences import fetch_species

from ttf.relational_external import shard_items

FIELDS = ["species", "source_key", "latitude", "longitude", "priority_rank", "priority_sha256"]


def _worker(species: str, out: mp.Queue) -> None:
    try:
        retained, ledger = fetch_species(species)
        out.put(("ok", retained, ledger))
    except BaseException as exc:
        out.put((
            "error",
            [],
            {
                "species": species,
                "status": "REQUEST_ERROR",
                "error": f"{type(exc).__name__}: {exc}",
            },
        ))


def bounded_fetch(species: str, timeout_seconds: int) -> tuple[list[dict[str, object]], dict[str, object]]:
    ctx = mp.get_context("spawn")
    out: mp.Queue = ctx.Queue(maxsize=1)
    process = ctx.Process(target=_worker, args=(species, out), daemon=True)
    process.start()
    process.join(timeout=float(timeout_seconds))
    if process.is_alive():
        process.terminate()
        process.join(timeout=15)
        if process.is_alive():
            process.kill()
            process.join(timeout=15)
        return [], {
            "species": species,
            "status": "REQUEST_ERROR",
            "error": f"PER_SPECIES_TRANSPORT_TIMEOUT_AFTER_{int(timeout_seconds)}S",
            "transport_timeout": True,
        }
    try:
        _, retained, ledger = out.get(timeout=10)
    except queue.Empty:
        return [], {
            "species": species,
            "status": "REQUEST_ERROR",
            "error": f"WORKER_EXIT_WITHOUT_RESULT_CODE_{process.exitcode}",
            "transport_timeout": False,
        }
    return list(retained), dict(ledger)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--candidates", type=Path, required=True)
    ap.add_argument("--output-csv", type=Path, required=True)
    ap.add_argument("--output-ledger", type=Path, required=True)
    ap.add_argument("--batch-index", type=int, required=True)
    ap.add_argument("--batches", type=int, default=250)
    ap.add_argument("--per-species-timeout-seconds", type=int, default=2700)
    args = ap.parse_args()

    rows = list(csv.DictReader(args.candidates.open(encoding="utf-8")))
    names = [str(row["species"]).strip() for row in rows]
    if len(names) != 1000 or len(set(names)) != 1000:
        raise RuntimeError("expected exact frozen Study-C 1000-species candidate set")
    if args.batches != 250:
        raise RuntimeError("bounded Study-C execution requires exactly 250 batches")
    if args.per_species_timeout_seconds != 2700:
        raise RuntimeError("bounded Study-C execution requires exact 2700-second per-species timeout")

    batch = list(shard_items(names, shard_index=args.batch_index, shards=args.batches))
    if len(batch) != 4:
        raise RuntimeError(f"expected exactly four species in batch {args.batch_index}, found {len(batch)}")

    retained_rows: list[dict[str, object]] = []
    ledgers: list[dict[str, object]] = []
    for ordinal, species in enumerate(batch):
        retained, ledger = bounded_fetch(species, args.per_species_timeout_seconds)
        retained_rows.extend(retained)
        ledger = dict(ledger)
        ledger["bounded_transport_batch_index"] = int(args.batch_index)
        ledger["bounded_transport_batch_ordinal"] = int(ordinal)
        ledger["per_species_wall_timeout_seconds"] = int(args.per_species_timeout_seconds)
        ledgers.append(ledger)
        print(json.dumps({
            "batch_index": args.batch_index,
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
        "schema": "ttf_relational_historical_occurrence_bounded_batch_v0.2",
        "batch_index": int(args.batch_index),
        "batches": int(args.batches),
        "candidate_species_total": 1000,
        "species_requested": batch,
        "retained_occurrence_rows": len(retained_rows),
        "species": ledgers,
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
