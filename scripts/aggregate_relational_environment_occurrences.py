#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from collections import Counter
from pathlib import Path


def sha256_path(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--input-dir", type=Path, required=True)
    ap.add_argument("--output-csv", type=Path, required=True)
    ap.add_argument("--output-ledger", type=Path, required=True)
    args = ap.parse_args()

    csv_paths = sorted(args.input_dir.glob("occurrences-*.csv"))
    ledger_paths = sorted(args.input_dir.glob("ledger-*.json"))
    if not csv_paths or not ledger_paths:
        raise RuntimeError("missing environmental occurrence shards")

    rows = []
    for path in csv_paths:
        rows.extend(csv.DictReader(path.open(encoding="utf-8")))
    species_ledgers = []
    for path in ledger_paths:
        payload = json.loads(path.read_text())
        if payload.get("schema") != "ttf_relational_environment_occurrence_shard_v0.2":
            raise RuntimeError(f"unexpected shard schema: {path}")
        if any(bool(v) for v in payload["response_firewall"].values()):
            raise RuntimeError("response firewall is open")
        species_ledgers.extend(payload["species"])

    names = [str(row["species"]) for row in species_ledgers]
    if len(names) != 1000 or len(set(names)) != 1000:
        raise RuntimeError("environment acquisition did not cover exact 1000 species")
    row_keys = [(row["species"], int(row["source_key"])) for row in rows]
    if len(row_keys) != len(set(row_keys)):
        raise RuntimeError("duplicate species/source_key rows across shards")
    rows.sort(key=lambda row: (row["species"], int(row["priority_rank"]), int(row["source_key"])))

    args.output_csv.parent.mkdir(parents=True, exist_ok=True)
    fields = ["species", "source_key", "latitude", "longitude", "priority_rank", "priority_sha256"]
    with args.output_csv.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)

    counts = Counter(str(row["status"]) for row in species_ledgers)
    passed = sum(1 for row in species_ledgers if row["status"] == "PASS_OCCURRENCE_GEOMETRY")
    payload = {
        "schema": "ttf_relational_environment_occurrence_acquisition_v0.2",
        "status": "PASS_TO_CHELSA_EXTRACTION" if passed >= 500 else "NOT_EVALUABLE_OCCURRENCE_GEOMETRY",
        "species": 1000,
        "species_passing_ge_30": passed,
        "status_counts": dict(sorted(counts.items())),
        "retained_rows": len(rows),
        "occurrence_csv_sha256": sha256_path(args.output_csv),
        "response_firewall": {
            "study_B_sequence_identity_opened": False,
            "study_B_pairwise_genetic_distances_opened": False,
            "study_B_T_st_computed": False,
            "study_B_beta_R_computed": False,
        },
    }
    args.output_ledger.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    print(json.dumps({k: payload[k] for k in ("status", "species_passing_ge_30", "status_counts", "retained_rows")}, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
