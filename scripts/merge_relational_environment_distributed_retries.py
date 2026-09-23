#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from collections import Counter
from pathlib import Path


FIELDS = ["species", "source_key", "latitude", "longitude", "priority_rank", "priority_sha256"]


def digest_species(names: list[str]) -> str:
    return hashlib.sha256(("\n".join(sorted(names)) + "\n").encode()).hexdigest()


def load_base(input_dir: Path):
    species: dict[str, dict] = {}
    rows_by_species: dict[str, list[dict[str, str]]] = {}
    for path in sorted(input_dir.glob("ledger-*.json")):
        payload = json.loads(path.read_text())
        if payload.get("schema") != "ttf_relational_environment_occurrence_shard_v0.2":
            raise RuntimeError(f"unexpected base shard schema: {path}")
        if any(bool(v) for v in payload["response_firewall"].values()):
            raise RuntimeError("base response firewall is open")
        for row in payload["species"]:
            name = str(row["species"])
            if name in species:
                raise RuntimeError(f"duplicate base species: {name}")
            species[name] = dict(row)
    for path in sorted(input_dir.glob("occurrences-*.csv")):
        for row in csv.DictReader(path.open(encoding="utf-8")):
            rows_by_species.setdefault(str(row["species"]), []).append(dict(row))
    if len(species) != 1000:
        raise RuntimeError(f"expected exact 1000 base species, found {len(species)}")
    return species, rows_by_species


def load_distributed_retry(retry_dir: Path):
    ledgers = sorted(retry_dir.glob("ledger-retry-*.json"))
    csvs = {p.stem.removeprefix("occurrences-retry-"): p for p in retry_dir.glob("occurrences-retry-*.csv")}
    if len(ledgers) != 235:
        raise RuntimeError(f"expected 235 distributed retry ledgers, found {len(ledgers)}")

    species_rows: dict[str, dict] = {}
    occurrence_rows: dict[str, list[dict[str, str]]] = {}
    attempt_history: dict[str, list[str]] = {}
    seen_batches = set()

    for path in ledgers:
        payload = json.loads(path.read_text())
        if payload.get("schema") != "ttf_relational_environment_distributed_retry_batch_v0.1":
            raise RuntimeError(f"unexpected distributed retry schema: {path}")
        if any(bool(v) for v in payload["response_firewall"].values()):
            raise RuntimeError("distributed retry firewall is open")
        batch_index = int(payload["batch_index"])
        if batch_index in seen_batches:
            raise RuntimeError(f"duplicate distributed retry batch: {batch_index}")
        seen_batches.add(batch_index)
        csv_path = csvs.get(str(batch_index))
        if csv_path is None:
            raise RuntimeError(f"missing occurrence CSV for distributed retry batch {batch_index}")
        batch_occurrences: dict[str, list[dict[str, str]]] = {}
        for row in csv.DictReader(csv_path.open(encoding="utf-8")):
            batch_occurrences.setdefault(str(row["species"]), []).append(dict(row))

        for row in payload["species"]:
            name = str(row["species"])
            if name in species_rows:
                raise RuntimeError(f"duplicate retry species: {name}")
            species_rows[name] = dict(row)
            occurrence_rows[name] = batch_occurrences.get(name, [])
        for name, hist in payload.get("attempt_history", {}).items():
            if name in attempt_history:
                raise RuntimeError(f"duplicate retry history species: {name}")
            attempt_history[str(name)] = list(map(str, hist))

    if seen_batches != set(range(235)):
        missing = sorted(set(range(235)) - seen_batches)
        extra = sorted(seen_batches - set(range(235)))
        raise RuntimeError(f"distributed retry batch coverage drift: missing={missing}, extra={extra}")
    return species_rows, occurrence_rows, attempt_history


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--base-dir", type=Path, required=True)
    ap.add_argument("--retry-dir", type=Path, required=True)
    ap.add_argument("--preretry-audit", type=Path, required=True)
    ap.add_argument("--fallback-trigger", type=Path, required=True)
    ap.add_argument("--output-csv", type=Path, required=True)
    ap.add_argument("--output-ledger", type=Path, required=True)
    ap.add_argument("--output-retry-receipt", type=Path, required=True)
    args = ap.parse_args()

    audit = json.loads(args.preretry_audit.read_text())
    trigger = json.loads(args.fallback_trigger.read_text())
    if audit.get("schema") != "ttf_relational_environment_transport_preretry_audit_v0.7":
        raise RuntimeError("unexpected pre-retry audit schema")
    if audit.get("status") != "RESPONSE_BLIND_PRE_RETRY_STATUS_AUDIT":
        raise RuntimeError("pre-retry audit status drift")
    if trigger.get("schema") != "ttf_relational_environment_request_error_fallback_v0.1":
        raise RuntimeError("unexpected distributed fallback trigger schema")
    if trigger.get("status") != "FROZEN_DISTRIBUTED_TECHNICAL_RETRY_BEFORE_RELATION_RESULT":
        raise RuntimeError("distributed fallback trigger is not frozen")
    if any(bool(v) for v in trigger["response_firewall"].values()):
        raise RuntimeError("distributed fallback trigger firewall is open")

    frozen_request_errors = list(map(str, audit["request_error_species"]))
    if len(frozen_request_errors) != 704 or len(set(frozen_request_errors)) != 704:
        raise RuntimeError("pre-retry request-error species count drift")
    frozen_digest = digest_species(frozen_request_errors)
    if frozen_digest != trigger["request_error_species_sorted_sha256"]:
        raise RuntimeError("pre-retry request-error species digest drift")

    base_species, rows_by_species = load_base(args.base_dir)
    base_request_errors = sorted(name for name, row in base_species.items() if row["status"] == "REQUEST_ERROR")
    if base_request_errors != sorted(frozen_request_errors):
        raise RuntimeError("base REQUEST_ERROR species set does not match frozen pre-retry audit")

    retry_species, retry_rows, attempt_history = load_distributed_retry(args.retry_dir)
    if sorted(retry_species) != sorted(frozen_request_errors):
        raise RuntimeError("distributed retry species set does not match frozen 704-species audit")

    for name in frozen_request_errors:
        if base_species[name]["status"] != "REQUEST_ERROR":
            raise RuntimeError(f"attempt to replace non-REQUEST_ERROR base species: {name}")
        base_species[name] = retry_species[name]
        rows_by_species[name] = retry_rows.get(name, [])

    unresolved = sorted(name for name, row in base_species.items() if row["status"] == "REQUEST_ERROR")
    all_rows = []
    for name in sorted(rows_by_species):
        all_rows.extend(rows_by_species[name])
    all_rows.sort(key=lambda row: (row["species"], int(row["priority_rank"]), int(row["source_key"])))

    args.output_csv.parent.mkdir(parents=True, exist_ok=True)
    with args.output_csv.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDS)
        writer.writeheader()
        writer.writerows(all_rows)

    counts = Counter(str(row["status"]) for row in base_species.values())
    passed = counts.get("PASS_OCCURRENCE_GEOMETRY", 0)
    status = (
        "INCOMPLETE_TECHNICAL_EXECUTION"
        if unresolved else
        ("PASS_TO_CHELSA_EXTRACTION" if passed >= 500 else "NOT_EVALUABLE_OCCURRENCE_GEOMETRY")
    )
    ledger = {
        "schema": "ttf_relational_environment_occurrence_acquisition_v0.2",
        "status": status,
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
    args.output_ledger.write_text(json.dumps(ledger, indent=2, sort_keys=True) + "\n")

    attempt_counts = Counter(len(hist) for hist in attempt_history.values())
    retry_receipt = {
        "schema": "ttf_relational_environment_distributed_retry_merge_v0.1",
        "status": (
            "PASS_ZERO_REQUEST_ERROR_DISTRIBUTED_MERGE"
            if not unresolved else
            "INCOMPLETE_DISTRIBUTED_RETRY_MERGE"
        ),
        "pre_retry_request_error_species": 704,
        "pre_retry_species_digest_sha256": frozen_digest,
        "final_request_error_count": len(unresolved),
        "final_status_counts": dict(sorted(counts.items())),
        "attempt_length_counts": {str(k): v for k, v in sorted(attempt_counts.items())},
        "response_firewall": {
            "Study_B_sequence_identity_opened": False,
            "Study_B_pairwise_genetic_distances_opened": False,
            "Study_B_T_st_computed": False,
            "Study_B_beta_R_computed": False,
        },
    }
    args.output_retry_receipt.write_text(json.dumps(retry_receipt, indent=2, sort_keys=True) + "\n")
    print(json.dumps({
        "status": status,
        "retry_status": retry_receipt["status"],
        "species_passing_ge_30": passed,
        "request_errors": len(unresolved),
        "status_counts": dict(sorted(counts.items())),
    }, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
