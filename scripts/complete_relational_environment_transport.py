#!/usr/bin/env python3
from __future__ import annotations

import argparse
import concurrent.futures
import csv
import json
from collections import Counter
from pathlib import Path

try:
    from scripts.acquire_relational_environment_occurrences import fetch_species
except ModuleNotFoundError:
    from acquire_relational_environment_occurrences import fetch_species


FIELDS = ["species", "source_key", "latitude", "longitude", "priority_rank", "priority_sha256"]
RETRY_WORKERS = 20


def load_base(input_dir: Path):
    species_ledger = {}
    rows_by_species = {}
    for path in sorted(input_dir.glob("ledger-*.json")):
        payload = json.loads(path.read_text())
        if payload.get("schema") != "ttf_relational_environment_occurrence_shard_v0.2":
            raise RuntimeError(f"unexpected shard schema: {path}")
        if any(bool(v) for v in payload["response_firewall"].values()):
            raise RuntimeError("base response firewall is open")
        for row in payload["species"]:
            name = str(row["species"])
            if name in species_ledger:
                raise RuntimeError(f"duplicate base species: {name}")
            species_ledger[name] = dict(row)
    for path in sorted(input_dir.glob("occurrences-*.csv")):
        for row in csv.DictReader(path.open(encoding="utf-8")):
            rows_by_species.setdefault(str(row["species"]), []).append(dict(row))
    if len(species_ledger) != 1000:
        raise RuntimeError(f"expected exact 1000 base species, found {len(species_ledger)}")
    return species_ledger, rows_by_species


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--input-dir", type=Path, required=True)
    ap.add_argument("--max-rounds", type=int, default=3)
    ap.add_argument("--output-csv", type=Path, required=True)
    ap.add_argument("--output-ledger", type=Path, required=True)
    ap.add_argument("--output-retry-receipt", type=Path, required=True)
    args = ap.parse_args()
    if args.max_rounds != 3:
        raise RuntimeError("frozen transport retry contract requires exactly 3 maximum rounds")

    species, rows_by_species = load_base(args.input_dir)
    history = []
    for round_id in range(1, args.max_rounds + 1):
        pending = sorted(name for name, row in species.items() if row["status"] == "REQUEST_ERROR")
        history.append({"round": round_id, "request_errors_before": len(pending)})
        if not pending:
            history[-1]["request_errors_after"] = 0
            break
        def retry_one(name: str):
            try:
                retained, ledger = fetch_species(name)
            except Exception as exc:
                retained, ledger = [], {
                    "species": name,
                    "status": "REQUEST_ERROR",
                    "error": f"{type(exc).__name__}: {exc}",
                }
            return name, retained, dict(ledger)

        completed: dict[str, tuple[list[dict[str, object]], dict[str, object]]] = {}
        with concurrent.futures.ThreadPoolExecutor(max_workers=RETRY_WORKERS) as pool:
            future_by_name = {pool.submit(retry_one, name): name for name in pending}
            for i, future in enumerate(concurrent.futures.as_completed(future_by_name), start=1):
                name, retained, ledger = future.result()
                completed[name] = (retained, ledger)
                if i % 20 == 0 or i == len(pending):
                    print(json.dumps({
                        "retry_round": round_id,
                        "completed": i,
                        "total": len(pending),
                        "workers": RETRY_WORKERS,
                    }))

        for name in pending:
            retained, ledger = completed[name]
            species[name] = ledger
            rows_by_species[name] = retained
        after = sum(row["status"] == "REQUEST_ERROR" for row in species.values())
        history[-1]["request_errors_after"] = after
        if after == 0:
            break

    unresolved = sorted(name for name, row in species.items() if row["status"] == "REQUEST_ERROR")
    all_rows = []
    for name in sorted(rows_by_species):
        all_rows.extend(rows_by_species[name])
    all_rows.sort(key=lambda row: (row["species"], int(row["priority_rank"]), int(row["source_key"])))

    args.output_csv.parent.mkdir(parents=True, exist_ok=True)
    with args.output_csv.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDS)
        writer.writeheader()
        writer.writerows(all_rows)

    counts = Counter(str(row["status"]) for row in species.values())
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

    receipt = {
        "schema": "ttf_relational_environment_transport_retry_receipt_v0.1",
        "status": (
            "PASS_ZERO_REQUEST_ERROR_TRANSPORT"
            if not unresolved else
            "INCOMPLETE_TECHNICAL_EXECUTION_AFTER_FROZEN_RETRIES"
        ),
        "maximum_retry_rounds": 3,
        "retry_workers": RETRY_WORKERS,
        "retry_history": history,
        "final_request_error_count": len(unresolved),
        "final_status_counts": dict(sorted(counts.items())),
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
    }, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
