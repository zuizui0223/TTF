#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

FIELDS = ["species", "source_key", "latitude", "longitude", "priority_rank", "priority_sha256"]


def load_occurrence_rows(path: Path) -> dict[str, list[dict[str, str]]]:
    out: dict[str, list[dict[str, str]]] = {}
    if not path.is_file():
        return out
    for row in csv.DictReader(path.open(encoding="utf-8")):
        out.setdefault(str(row["species"]), []).append(dict(row))
    return out


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--trigger", type=Path, required=True)
    ap.add_argument("--repair-dir", type=Path, required=True)
    ap.add_argument("--recovery-csv", type=Path, required=True)
    ap.add_argument("--recovery-ledger", type=Path, required=True)
    ap.add_argument("--output-csv", type=Path, required=True)
    ap.add_argument("--output-ledger", type=Path, required=True)
    args = ap.parse_args()

    trigger = json.loads(args.trigger.read_text())
    if trigger.get("schema") != "ttf_relational_environment_cutover_timeout_repair_v0.1":
        raise RuntimeError("unexpected frozen cutover-repair trigger")
    entries = sorted(trigger["entries"], key=lambda x: int(x["repair_index"]))
    if len(entries) != 8:
        raise RuntimeError("expected exact eight frozen cutover-repair entries")
    expected_species = [str(x["species"]) for x in entries]
    expected_indices = [int(x["candidate_index"]) for x in entries]

    original_ledger: dict[str, dict] = {}
    original_rows: dict[str, list[dict[str, str]]] = {}
    for path in sorted(args.repair_dir.glob("ledger-cutover-repair-*.json")):
        payload = json.loads(path.read_text())
        if payload.get("schema") != "ttf_relational_environment_occurrence_shard_v0.2":
            raise RuntimeError(f"unexpected original repair schema: {path}")
        if any(bool(v) for v in payload["response_firewall"].values()):
            raise RuntimeError("original repair response firewall is open")
        for row in payload["species"]:
            name = str(row["species"])
            if name in original_ledger:
                raise RuntimeError(f"duplicate original repair species: {name}")
            original_ledger[name] = dict(row)
    for path in sorted(args.repair_dir.glob("occurrences-cutover-repair-*.csv")):
        for name, rows in load_occurrence_rows(path).items():
            original_rows.setdefault(name, []).extend(rows)

    if set(original_ledger) != set(expected_species):
        raise RuntimeError(
            f"original repair species drift: expected={expected_species}, got={sorted(original_ledger)}"
        )

    recovery = json.loads(args.recovery_ledger.read_text())
    if recovery.get("schema") != "ttf_relational_environment_rate_limit_recovery_result_v0.1":
        raise RuntimeError("unexpected rate-limit recovery result schema")
    if any(bool(v) for v in recovery["response_firewall"].values()):
        raise RuntimeError("rate-limit recovery response firewall is open")
    if int(recovery["request_error_count"]) != 0:
        raise RuntimeError("rate-limit recovery still contains REQUEST_ERROR")
    recovery_ledger = {str(row["species"]): dict(row) for row in recovery["species"]}
    frozen_recovery = {
        "Pyrrhosoma nymphula",
        "Lithobates clamitans",
        "Diarsia rubi",
        "Mythimna impura",
    }
    if set(recovery_ledger) != frozen_recovery:
        raise RuntimeError("rate-limit recovery species drift")
    recovery_rows = load_occurrence_rows(args.recovery_csv)

    final_ledger = {}
    final_rows = {}
    for name in expected_species:
        if name in recovery_ledger:
            final_ledger[name] = recovery_ledger[name]
            final_rows[name] = recovery_rows.get(name, [])
        else:
            if original_ledger[name]["status"] == "REQUEST_ERROR":
                raise RuntimeError(f"unrecovered original REQUEST_ERROR: {name}")
            final_ledger[name] = original_ledger[name]
            final_rows[name] = original_rows.get(name, [])

    if any(row["status"] == "REQUEST_ERROR" for row in final_ledger.values()):
        raise RuntimeError("final eight-species repair contains REQUEST_ERROR")

    all_rows = []
    for name in expected_species:
        all_rows.extend(final_rows.get(name, []))
    all_rows.sort(
        key=lambda row: (
            expected_species.index(str(row["species"])),
            int(row["priority_rank"]),
            int(row["source_key"]),
        )
    )
    args.output_csv.parent.mkdir(parents=True, exist_ok=True)
    with args.output_csv.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDS)
        writer.writeheader()
        writer.writerows(all_rows)

    payload = {
        "schema": "ttf_relational_environment_occurrence_shard_v0.2",
        "shard_index": 0,
        "shards": 1,
        "species_requested": 8,
        "candidate_indices": expected_indices,
        "retained_occurrence_rows": len(all_rows),
        "species": [final_ledger[name] for name in expected_species],
        "final_cutover_timeout_repair": True,
        "response_firewall": {
            "study_B_sequence_identity_opened": False,
            "study_B_pairwise_genetic_distances_opened": False,
            "study_B_T_st_computed": False,
            "study_B_beta_R_computed": False,
        },
    }
    args.output_ledger.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    print(json.dumps({
        "species": 8,
        "request_errors": sum(x["status"] == "REQUEST_ERROR" for x in payload["species"]),
        "statuses": {x["species"]: x["status"] for x in payload["species"]},
    }, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
