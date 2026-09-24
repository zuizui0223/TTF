#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import json
import multiprocessing as mp
from pathlib import Path
from queue import Empty

try:
    from scripts.acquire_relational_environment_occurrences import fetch_species
except ModuleNotFoundError:
    from acquire_relational_environment_occurrences import fetch_species


SCHEMA = "ttf_relational_historical_occurrence_shard_v0.1"
CANDIDATE_TOTAL = 1000
FROZEN_BATCH_SIZE = 4
FROZEN_BATCH_COUNT = 250
FROZEN_SPECIES_TIMEOUT_SECONDS = 3300


def batch_species(names: list[str], batch_index: int, batch_size: int = FROZEN_BATCH_SIZE) -> list[str]:
    if len(names) != CANDIDATE_TOTAL or len(set(names)) != CANDIDATE_TOTAL:
        raise RuntimeError("Study C repair requires exact frozen 1000-species candidate set")
    if batch_size != FROZEN_BATCH_SIZE:
        raise RuntimeError("Study C repair batch size drift")
    if not 0 <= int(batch_index) < FROZEN_BATCH_COUNT:
        raise ValueError("invalid Study C repair batch index")
    start = int(batch_index) * batch_size
    return list(names[start : start + batch_size])


def timeout_ledger(species: str, timeout_seconds: int) -> dict[str, object]:
    return {
        "species": str(species),
        "status": "REQUEST_ERROR",
        "error": f"TRANSPORT_HARD_TIMEOUT_AFTER_{int(timeout_seconds)}_SECONDS",
    }


def _worker(species: str, queue: mp.Queue) -> None:
    try:
        retained, ledger = fetch_species(species)
    except Exception as exc:
        retained, ledger = [], {
            "species": species,
            "status": "REQUEST_ERROR",
            "error": f"{type(exc).__name__}: {exc}",
        }
    queue.put((retained, ledger))


def fetch_with_hard_timeout(species: str, timeout_seconds: int) -> tuple[list[dict[str, object]], dict[str, object]]:
    if int(timeout_seconds) != FROZEN_SPECIES_TIMEOUT_SECONDS:
        raise RuntimeError("Study C per-species timeout drift")
    ctx = mp.get_context("fork")
    queue = ctx.Queue(maxsize=1)
    proc = ctx.Process(target=_worker, args=(species, queue))
    proc.start()
    proc.join(timeout=float(timeout_seconds))
    if proc.is_alive():
        proc.terminate()
        proc.join(timeout=10)
        return [], timeout_ledger(species, timeout_seconds)
    try:
        retained, ledger = queue.get(timeout=5)
    except Empty:
        return [], {
            "species": species,
            "status": "REQUEST_ERROR",
            "error": f"TRANSPORT_CHILD_EXIT_WITHOUT_RESULT_EXITCODE_{proc.exitcode}",
        }
    return list(retained), dict(ledger)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--candidates", type=Path, required=True)
    ap.add_argument("--output-csv", type=Path, required=True)
    ap.add_argument("--output-ledger", type=Path, required=True)
    ap.add_argument("--batch-index", type=int, required=True)
    ap.add_argument("--batch-size", type=int, default=FROZEN_BATCH_SIZE)
    ap.add_argument(
        "--species-timeout-seconds",
        type=int,
        default=FROZEN_SPECIES_TIMEOUT_SECONDS,
    )
    args = ap.parse_args()

    rows = list(csv.DictReader(args.candidates.open(encoding="utf-8")))
    names = [str(row["species"]).strip() for row in rows]
    selected = batch_species(names, args.batch_index, args.batch_size)

    retained_rows: list[dict[str, object]] = []
    ledgers: list[dict[str, object]] = []
    for species in selected:
        retained, ledger = fetch_with_hard_timeout(
            species,
            int(args.species_timeout_seconds),
        )
        retained_rows.extend(retained)
        ledgers.append(ledger)
        print(json.dumps({
            "batch_index": int(args.batch_index),
            "species": species,
            "status": ledger["status"],
        }, sort_keys=True), flush=True)

    args.output_csv.parent.mkdir(parents=True, exist_ok=True)
    fields = [
        "species",
        "source_key",
        "latitude",
        "longitude",
        "priority_rank",
        "priority_sha256",
    ]
    with args.output_csv.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(retained_rows)

    payload = {
        "schema": SCHEMA,
        "execution_version": "v0.2-singleton-safe",
        "batch_index": int(args.batch_index),
        "shard_index": int(args.batch_index),
        "shards": FROZEN_BATCH_COUNT,
        "candidate_species_total": CANDIDATE_TOTAL,
        "species_requested": len(selected),
        "species_names": selected,
        "retained_occurrence_rows": len(retained_rows),
        "species": ledgers,
        "per_species_hard_timeout_seconds": FROZEN_SPECIES_TIMEOUT_SECONDS,
        "inheritance": (
            "Exact Study-B GBIF fetch_species query/selection logic reused unchanged; "
            "outer hard timeout only guarantees a terminal technical ledger."
        ),
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
