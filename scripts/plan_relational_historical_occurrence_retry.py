#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
from collections import Counter
from pathlib import Path


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--input-dir", type=Path, required=True)
    ap.add_argument("--candidates", type=Path, required=True)
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()

    import csv
    candidate_rows = list(csv.DictReader(args.candidates.open(encoding="utf-8")))
    names = [str(row["species"]).strip() for row in candidate_rows]
    if len(names) != 1000 or len(set(names)) != 1000:
        raise RuntimeError("expected exact frozen Study-C 1000-species candidate set")

    ledgers = sorted(args.input_dir.glob("ledger-bounded-*.json"))
    if len(ledgers) != 250:
        raise RuntimeError(f"expected exact 250 bounded Study-C ledgers, found {len(ledgers)}")

    by_name: dict[str, dict] = {}
    seen_batches: set[int] = set()
    for path in ledgers:
        p = json.loads(path.read_text())
        if p.get("schema") != "ttf_relational_historical_occurrence_bounded_batch_v0.2":
            raise RuntimeError(f"unexpected bounded Study-C ledger: {path}")
        if p.get("scientific_query_change") is not False:
            raise RuntimeError("bounded Study-C acquisition changed scientific query")
        if any(bool(v) for v in p["response_firewall"].values()):
            raise RuntimeError("bounded Study-C response firewall is open")
        batch = int(p["batch_index"])
        if batch in seen_batches:
            raise RuntimeError(f"duplicate bounded Study-C batch {batch}")
        seen_batches.add(batch)
        expected = [name for i, name in enumerate(names) if i % 250 == batch]
        observed = list(map(str, p["species_requested"]))
        if observed != expected:
            raise RuntimeError(f"bounded Study-C batch assignment drift at {batch}")
        for row in p["species"]:
            name = str(row["species"])
            if name in by_name:
                raise RuntimeError(f"duplicate bounded Study-C species {name}")
            by_name[name] = dict(row)

    if seen_batches != set(range(250)):
        raise RuntimeError("bounded Study-C batch coverage drift")
    if set(by_name) != set(names):
        raise RuntimeError("bounded Study-C species universe drift")

    request_errors = [name for name in names if by_name[name]["status"] == "REQUEST_ERROR"]
    counts = Counter(str(by_name[name]["status"]) for name in names)
    payload = {
        "schema": "ttf_relational_historical_occurrence_retry_plan_v0.2",
        "status": "RETRY_EXACT_INITIAL_REQUEST_ERRORS" if request_errors else "NO_RETRY_NEEDED",
        "candidate_species": 1000,
        "initial_batches": 250,
        "request_error_count": len(request_errors),
        "request_error_species": request_errors,
        "status_counts": dict(sorted(counts.items())),
        "retry_batch_size": 4,
        "retry_batch_count": (len(request_errors) + 3) // 4,
        "maximum_retry_rounds": 1,
        "scientific_query_change": False,
        "relation_result_seen": False,
        "genetic_response_used": False,
        "response_firewall": {
            "Study_C_sequence_identity_opened": False,
            "Study_C_pairwise_genetic_distances_opened": False,
            "Study_C_T_st_computed": False,
            "Study_C_beta_hist_computed": False,
        },
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    print(json.dumps({
        "status": payload["status"],
        "request_error_count": payload["request_error_count"],
        "retry_batch_count": payload["retry_batch_count"],
        "status_counts": payload["status_counts"],
    }, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
