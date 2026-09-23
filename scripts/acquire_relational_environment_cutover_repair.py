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


FIELDS = ["species", "source_key", "latitude", "longitude", "priority_rank", "priority_sha256"]


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--trigger", type=Path, required=True)
    ap.add_argument("--repair-index", type=int, required=True)
    ap.add_argument("--output-csv", type=Path, required=True)
    ap.add_argument("--output-ledger", type=Path, required=True)
    args = ap.parse_args()

    trigger = json.loads(args.trigger.read_text())
    allowed = {
        "ttf_relational_environment_cutover_timeout_repair_v0.1":
            "FROZEN_CUTOVER_TIMEOUT_REPAIR_BEFORE_RELATION_RESULT",
        "ttf_relational_environment_cutover_timeout_repair_long_v0.1":
            "FROZEN_LONG_TIMEOUT_RETRY_AFTER_SINGLETON_TIMEOUT",
    }
    schema = trigger.get("schema")
    if schema not in allowed:
        raise RuntimeError("unexpected Study-B cutover repair trigger schema")
    if trigger.get("status") != allowed[schema]:
        raise RuntimeError("Study-B cutover repair trigger is not frozen")
    if trigger.get("query_or_scientific_contract_change") is not False:
        raise RuntimeError("cutover repair changed scientific query contract")
    if any(bool(v) for v in trigger["response_firewall"].values()):
        raise RuntimeError("cutover repair response firewall is open")

    entries = list(trigger["entries"])
    if not 0 <= args.repair_index < len(entries):
        raise RuntimeError("repair index outside frozen trigger")
    entry = entries[args.repair_index]
    if int(entry["repair_index"]) != args.repair_index:
        raise RuntimeError("repair index drift")
    species = str(entry["species"])
    candidate_index = int(entry["candidate_index"])

    try:
        retained, ledger = fetch_species(species)
    except Exception as exc:
        retained, ledger = [], {
            "species": species,
            "status": "REQUEST_ERROR",
            "error": f"{type(exc).__name__}: {exc}",
        }

    args.output_csv.parent.mkdir(parents=True, exist_ok=True)
    with args.output_csv.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDS)
        writer.writeheader()
        writer.writerows(retained)

    payload = {
        "schema": "ttf_relational_environment_occurrence_shard_v0.2",
        "shard_index": args.repair_index,
        "shards": len(entries),
        "species_requested": 1,
        "candidate_indices": [candidate_index],
        "replaces_cutover_batch": int(entry["cutover_batch_index"]),
        "retained_occurrence_rows": len(retained),
        "species": [ledger],
        "response_firewall": {
            "study_B_sequence_identity_opened": False,
            "study_B_pairwise_genetic_distances_opened": False,
            "study_B_T_st_computed": False,
            "study_B_beta_R_computed": False,
        },
    }
    args.output_ledger.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
