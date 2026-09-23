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


def closed_firewall(payload: dict) -> None:
    firewall = payload.get("response_firewall")
    if not isinstance(firewall, dict) or not firewall:
        raise RuntimeError("missing response firewall")
    if any(bool(v) for v in firewall.values()):
        raise RuntimeError("response firewall is open")


def load_base(input_dir: Path) -> tuple[dict[str, dict], dict[str, list[dict[str, str]]]]:
    species: dict[str, dict] = {}
    rows_by_species: dict[str, list[dict[str, str]]] = {}
    for path in sorted(input_dir.glob("ledger-*.json")):
        payload = json.loads(path.read_text())
        if payload.get("schema") != "ttf_relational_environment_occurrence_shard_v0.2":
            raise RuntimeError(f"unexpected base shard schema: {path}")
        closed_firewall(payload)
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


def load_round3(input_dir: Path) -> tuple[dict[str, dict], dict[str, list[dict[str, str]]], dict[str, list[str]]]:
    ledgers = sorted(input_dir.glob("ledger-retry-*.json"))
    csvs = {
        p.stem.removeprefix("occurrences-retry-"): p
        for p in input_dir.glob("occurrences-retry-*.csv")
    }
    if len(ledgers) != 235:
        raise RuntimeError(f"expected 235 round-3 ledgers, found {len(ledgers)}")
    species_rows: dict[str, dict] = {}
    occurrence_rows: dict[str, list[dict[str, str]]] = {}
    history: dict[str, list[str]] = {}
    batches: set[int] = set()
    for path in ledgers:
        payload = json.loads(path.read_text())
        if payload.get("schema") != "ttf_relational_environment_distributed_retry_batch_v0.1":
            raise RuntimeError(f"unexpected round-3 schema: {path}")
        closed_firewall(payload)
        batch = int(payload["batch_index"])
        if batch in batches:
            raise RuntimeError(f"duplicate round-3 batch: {batch}")
        batches.add(batch)
        csv_path = csvs.get(str(batch))
        if csv_path is None:
            raise RuntimeError(f"missing round-3 CSV for batch {batch}")
        by_name: dict[str, list[dict[str, str]]] = {}
        for row in csv.DictReader(csv_path.open(encoding="utf-8")):
            by_name.setdefault(str(row["species"]), []).append(dict(row))
        for row in payload["species"]:
            name = str(row["species"])
            if name in species_rows:
                raise RuntimeError(f"duplicate round-3 species: {name}")
            species_rows[name] = dict(row)
            occurrence_rows[name] = by_name.get(name, [])
        for name, attempts in payload.get("attempt_history", {}).items():
            if str(name) in history:
                raise RuntimeError(f"duplicate round-3 attempt history: {name}")
            history[str(name)] = list(map(str, attempts))
    if batches != set(range(235)):
        raise RuntimeError("round-3 batch coverage drift")
    return species_rows, occurrence_rows, history


def load_round4(input_dir: Path) -> tuple[dict[str, dict], dict[str, list[dict[str, str]]]]:
    ledgers = sorted(input_dir.glob("ledger-recovery-v02-*.json"))
    csvs = {
        p.stem.removeprefix("occurrences-recovery-v02-"): p
        for p in input_dir.glob("occurrences-recovery-v02-*.csv")
    }
    if len(ledgers) != 186:
        raise RuntimeError(f"expected 186 round-4 ledgers, found {len(ledgers)}")
    species_rows: dict[str, dict] = {}
    occurrence_rows: dict[str, list[dict[str, str]]] = {}
    batches: set[int] = set()
    for path in ledgers:
        payload = json.loads(path.read_text())
        if payload.get("schema") != "ttf_relational_environment_request_error_recovery_batch_v0.2":
            raise RuntimeError(f"unexpected round-4 schema: {path}")
        if int(payload.get("retry_round", -1)) != 4:
            raise RuntimeError(f"unexpected round-4 retry index: {path}")
        if payload.get("scientific_query_change") is not False:
            raise RuntimeError("round-4 scientific query changed")
        if payload.get("relation_result_seen") is not False or payload.get("genetic_response_used") is not False:
            raise RuntimeError("round-4 result was not response-blind")
        closed_firewall(payload)
        batch = int(payload["batch_index"])
        if batch in batches:
            raise RuntimeError(f"duplicate round-4 batch: {batch}")
        batches.add(batch)
        csv_path = csvs.get(str(batch))
        if csv_path is None:
            raise RuntimeError(f"missing round-4 CSV for batch {batch}")
        by_name: dict[str, list[dict[str, str]]] = {}
        for row in csv.DictReader(csv_path.open(encoding="utf-8")):
            by_name.setdefault(str(row["species"]), []).append(dict(row))
        for row in payload["species"]:
            name = str(row["species"])
            if name in species_rows:
                raise RuntimeError(f"duplicate round-4 species: {name}")
            species_rows[name] = dict(row)
            occurrence_rows[name] = by_name.get(name, [])
    if batches != set(range(186)):
        raise RuntimeError("round-4 batch coverage drift")
    return species_rows, occurrence_rows


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--base-dir", type=Path, required=True)
    ap.add_argument("--round3-dir", type=Path, required=True)
    ap.add_argument("--round4-dir", type=Path, required=True)
    ap.add_argument("--preretry-audit", type=Path, required=True)
    ap.add_argument("--round3-audit", type=Path, required=True)
    ap.add_argument("--round4-audit", type=Path, required=True)
    ap.add_argument("--round4-contract", type=Path, required=True)
    ap.add_argument("--output-csv", type=Path, required=True)
    ap.add_argument("--output-ledger", type=Path, required=True)
    ap.add_argument("--output-retry-receipt", type=Path, required=True)
    args = ap.parse_args()

    preretry = json.loads(args.preretry_audit.read_text())
    round3_audit = json.loads(args.round3_audit.read_text())
    round4_audit = json.loads(args.round4_audit.read_text())
    contract = json.loads(args.round4_contract.read_text())

    if preretry.get("schema") != "ttf_relational_environment_transport_preretry_audit_v0.7":
        raise RuntimeError("unexpected pre-retry audit schema")
    if preretry.get("status") != "RESPONSE_BLIND_PRE_RETRY_STATUS_AUDIT":
        raise RuntimeError("pre-retry audit status drift")
    closed_firewall(preretry)

    if round3_audit.get("schema") != "ttf_relational_environment_request_error_fallback_audit_v0.1":
        raise RuntimeError("unexpected round-3 audit schema")
    if round3_audit.get("status") != "INCOMPLETE_DISTRIBUTED_RETRY":
        raise RuntimeError("round-3 audit status drift")
    if int(round3_audit.get("species", -1)) != 704:
        raise RuntimeError("round-3 audit species count drift")
    if int(round3_audit["status_counts"].get("REQUEST_ERROR", -1)) != 556:
        raise RuntimeError("round-3 unresolved count drift")
    closed_firewall(round3_audit)

    if contract.get("schema") != "ttf_relational_environment_request_error_recovery_v0.2":
        raise RuntimeError("unexpected round-4 contract schema")
    if contract.get("status") != "FROZEN_FINAL_TRANSPORT_ROUND_BEFORE_RELATION_RESULT":
        raise RuntimeError("round-4 contract is not frozen")
    if int(contract.get("retry_round", -1)) != 4:
        raise RuntimeError("round-4 contract retry index drift")
    if int(contract.get("source_fallback_run_id", -1)) != int(round3_audit["fallback_workflow_run_id"]):
        raise RuntimeError("round-4 contract source fallback drift")
    if int(contract.get("source_unresolved_request_error_count", -1)) != 556:
        raise RuntimeError("round-4 contract unresolved count drift")
    closed_firewall(contract)

    if round4_audit.get("schema") != "ttf_relational_environment_request_error_recovery_audit_v0.2":
        raise RuntimeError("unexpected round-4 audit schema")
    if round4_audit.get("status") not in {
        "PASS_FINAL_ROUND_ZERO_REQUEST_ERROR",
        "INCOMPLETE_TECHNICAL_EXECUTION_AFTER_FINAL_ROUND",
    }:
        raise RuntimeError("round-4 audit status drift")
    if int(round4_audit.get("retry_round", -1)) != 4:
        raise RuntimeError("round-4 audit retry index drift")
    if int(round4_audit.get("recovery_workflow_run_id", -1)) <= 0:
        raise RuntimeError("round-4 audit lacks producer run binding")
    if round4_audit.get("no_fifth_round_authorized") is not True:
        raise RuntimeError("round-4 audit does not forbid a fifth round")
    closed_firewall(round4_audit)

    original = list(map(str, preretry["request_error_species"]))
    if len(original) != 704 or len(set(original)) != 704:
        raise RuntimeError("pre-retry REQUEST_ERROR universe drift")
    if digest_species(original) != contract["source_original_request_error_species_sorted_sha256"]:
        raise RuntimeError("pre-retry REQUEST_ERROR digest drift")
    if digest_species(original) != round3_audit["species_sorted_sha256"]:
        raise RuntimeError("round-3 audit species digest drift")

    frozen_round4 = list(map(str, round3_audit["unresolved_request_error_species"]))
    if len(frozen_round4) != 556 or len(set(frozen_round4)) != 556:
        raise RuntimeError("round-4 frozen species universe drift")

    base_species, rows_by_species = load_base(args.base_dir)
    base_request_errors = sorted(
        name for name, row in base_species.items() if row["status"] == "REQUEST_ERROR"
    )
    if base_request_errors != sorted(original):
        raise RuntimeError("base REQUEST_ERROR species do not match pre-retry audit")

    round3_species, round3_rows, round3_history = load_round3(args.round3_dir)
    if sorted(round3_species) != sorted(original):
        raise RuntimeError("round-3 species set does not match frozen 704")
    for name in original:
        if base_species[name]["status"] != "REQUEST_ERROR":
            raise RuntimeError(f"round-3 attempted to replace non-REQUEST_ERROR base species: {name}")
        base_species[name] = dict(round3_species[name])
        rows_by_species[name] = list(round3_rows.get(name, []))

    round3_unresolved = sorted(
        name for name, row in base_species.items() if row["status"] == "REQUEST_ERROR"
    )
    if round3_unresolved != sorted(frozen_round4):
        raise RuntimeError("round-3 merged unresolved set does not match frozen round-3 audit")
    retry3_counts = Counter(str(round3_species[name]["status"]) for name in original)
    if dict(sorted(retry3_counts.items())) != dict(sorted(round3_audit["status_counts"].items())):
        raise RuntimeError("round-3 merged status counts do not match frozen audit")

    round4_species, round4_rows = load_round4(args.round4_dir)
    if sorted(round4_species) != sorted(frozen_round4):
        raise RuntimeError("round-4 species set does not match exact frozen 556")
    for name in frozen_round4:
        if base_species[name]["status"] != "REQUEST_ERROR":
            raise RuntimeError(f"round-4 attempted to overwrite resolved round-3 species: {name}")
        base_species[name] = dict(round4_species[name])
        rows_by_species[name] = list(round4_rows.get(name, []))

    unresolved = sorted(
        name for name, row in base_species.items() if row["status"] == "REQUEST_ERROR"
    )
    if unresolved != sorted(map(str, round4_audit["unresolved_request_error_species"])):
        raise RuntimeError("final unresolved set does not match round-4 audit")
    if len(unresolved) != int(round4_audit["request_error_count_after_final_round"]):
        raise RuntimeError("final unresolved count does not match round-4 audit")
    if (not unresolved) != (round4_audit["status"] == "PASS_FINAL_ROUND_ZERO_REQUEST_ERROR"):
        raise RuntimeError("round-4 zero-error status inconsistent with merged transport")

    all_rows: list[dict[str, str]] = []
    for name in sorted(rows_by_species):
        all_rows.extend(rows_by_species[name])
    all_rows.sort(key=lambda row: (row["species"], int(row["priority_rank"]), int(row["source_key"])))

    args.output_csv.parent.mkdir(parents=True, exist_ok=True)
    with args.output_csv.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDS)
        writer.writeheader()
        writer.writerows(all_rows)

    counts = Counter(str(row["status"]) for row in base_species.values())
    passed = int(counts.get("PASS_OCCURRENCE_GEOMETRY", 0))
    if unresolved:
        status = "INCOMPLETE_TECHNICAL_EXECUTION_AFTER_FINAL_ROUND"
    elif passed >= 500:
        status = "PASS_TO_CHELSA_EXTRACTION"
    else:
        status = "NOT_EVALUABLE_OCCURRENCE_GEOMETRY"

    ledger = {
        "schema": "ttf_relational_environment_occurrence_acquisition_v0.3",
        "status": status,
        "species": 1000,
        "species_passing_ge_30": passed,
        "status_counts": dict(sorted(counts.items())),
        "retained_rows": len(all_rows),
        "unresolved_request_error_species": unresolved,
        "retry_rounds_completed": 4,
        "no_fifth_round_authorized": True,
        "response_firewall": {
            "Study_B_sequence_identity_opened": False,
            "Study_B_pairwise_genetic_distances_opened": False,
            "Study_B_T_st_computed": False,
            "Study_B_beta_R_computed": False,
        },
    }
    args.output_ledger.write_text(json.dumps(ledger, indent=2, sort_keys=True) + "\n")

    history_counts = Counter(len(hist) for hist in round3_history.values())
    receipt = {
        "schema": "ttf_relational_environment_final_transport_merge_v0.2",
        "status": (
            "PASS_ZERO_REQUEST_ERROR_FINAL_MERGE"
            if not unresolved
            else "INCOMPLETE_TECHNICAL_EXECUTION_AFTER_FINAL_ROUND"
        ),
        "pre_retry_request_error_species": 704,
        "round3_resolved_species": 148,
        "round3_unresolved_species": 556,
        "round4_species": 556,
        "final_request_error_count": len(unresolved),
        "final_status_counts": dict(sorted(counts.items())),
        "round3_attempt_length_counts": {str(k): v for k, v in sorted(history_counts.items())},
        "no_fifth_round_authorized": True,
        "response_firewall": {
            "Study_B_sequence_identity_opened": False,
            "Study_B_pairwise_genetic_distances_opened": False,
            "Study_B_T_st_computed": False,
            "Study_B_beta_R_computed": False,
        },
    }
    args.output_retry_receipt.write_text(json.dumps(receipt, indent=2, sort_keys=True) + "\n")
    print(json.dumps({
        "status": status,
        "retry_status": receipt["status"],
        "species_passing_ge_30": passed,
        "request_errors": len(unresolved),
        "status_counts": dict(sorted(counts.items())),
    }, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
