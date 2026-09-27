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


def request_error_species(input_dir: Path) -> list[str]:
    out = []
    seen = set()
    for path in sorted(input_dir.glob("ledger-*.json")):
        payload = json.loads(path.read_text())
        if payload.get("schema") != "ttf_relational_environment_occurrence_shard_v0.2":
            raise RuntimeError(f"unexpected shard schema: {path}")
        if any(bool(v) for v in payload["response_firewall"].values()):
            raise RuntimeError("response firewall is open")
        for row in payload["species"]:
            name = str(row["species"])
            if name in seen:
                raise RuntimeError(f"duplicate species across shard ledgers: {name}")
            seen.add(name)
            if row["status"] == "REQUEST_ERROR":
                out.append(name)
    if len(seen) != 1000:
        raise RuntimeError(f"expected exact 1000 species before retry, found {len(seen)}")
    return out


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--input-dir", type=Path, required=True)
    ap.add_argument("--round", type=int, required=True, choices=(1, 2, 3))
    ap.add_argument("--output-csv", type=Path, required=True)
    ap.add_argument("--output-ledger", type=Path, required=True)
    args = ap.parse_args()

    names = request_error_species(args.input_dir)
    retained_rows = []
    ledgers = []
    for i, species in enumerate(names, start=1):
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
            print(json.dumps({"retry_round": args.round, "completed": i, "total": len(names)}))

    args.output_csv.parent.mkdir(parents=True, exist_ok=True)
    fields = ["species", "source_key", "latitude", "longitude", "priority_rank", "priority_sha256"]
    with args.output_csv.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(retained_rows)

    unresolved = sum(row["status"] == "REQUEST_ERROR" for row in ledgers)
    payload = {
        "schema": "ttf_relational_environment_transport_retry_v0.1",
        "round": args.round,
        "requested_species": len(names),
        "unresolved_request_errors": unresolved,
        "species": ledgers,
        "response_firewall": {
            "Study_B_sequence_identity_opened": False,
            "Study_B_pairwise_genetic_distances_opened": False,
            "Study_B_T_st_computed": False,
            "Study_B_beta_R_computed": False,
        },
    }
    args.output_ledger.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    print(json.dumps({
        "round": args.round,
        "requested_species": len(names),
        "unresolved_request_errors": unresolved,
    }, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
