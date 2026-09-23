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
    ap.add_argument("--manifest", type=Path, required=True)
    ap.add_argument("--batch-index", type=int, required=True)
    ap.add_argument("--output-csv", type=Path, required=True)
    ap.add_argument("--output-ledger", type=Path, required=True)
    args = ap.parse_args()

    manifest = json.loads(args.manifest.read_text())
    if manifest.get("schema") != "ttf_relational_environment_transport_cutover_v0.1":
        raise RuntimeError("unexpected Study-B cutover manifest schema")
    if manifest.get("status") != "FROZEN_FINAL_TRANSPORT_CUTOVER_BEFORE_RELATION_RESULT":
        raise RuntimeError("Study-B cutover manifest is not frozen")
    if manifest.get("query_or_scientific_contract_change") is not False:
        raise RuntimeError("Study-B cutover changed the scientific query contract")
    if any(bool(v) for v in manifest["response_firewall"].values()):
        raise RuntimeError("Study-B cutover response firewall is open")

    batches = list(manifest["batches"])
    if not 0 <= args.batch_index < len(batches):
        raise RuntimeError("batch index outside frozen cutover manifest")
    batch = batches[args.batch_index]
    if int(batch["batch_index"]) != args.batch_index:
        raise RuntimeError("batch index drift in cutover manifest")
    names = list(map(str, batch["species"]))
    candidate_indices = list(map(int, batch["candidate_indices"]))
    if not (1 <= len(names) <= 4) or len(names) != len(candidate_indices):
        raise RuntimeError("invalid frozen cutover batch")
    if len(set(names)) != len(names) or len(set(candidate_indices)) != len(candidate_indices):
        raise RuntimeError("duplicate species/index inside frozen cutover batch")

    retained_rows = []
    ledgers = []
    for species in names:
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

    args.output_csv.parent.mkdir(parents=True, exist_ok=True)
    with args.output_csv.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDS)
        writer.writeheader()
        writer.writerows(retained_rows)

    payload = {
        "schema": "ttf_relational_environment_occurrence_shard_v0.2",
        "shard_index": args.batch_index,
        "shards": len(batches),
        "species_requested": len(names),
        "candidate_indices": candidate_indices,
        "retained_occurrence_rows": len(retained_rows),
        "species": ledgers,
        "cutover_manifest_schema": manifest["schema"],
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
