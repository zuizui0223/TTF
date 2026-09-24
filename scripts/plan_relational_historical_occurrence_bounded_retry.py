#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path


def main() -> int:
    ap=argparse.ArgumentParser()
    ap.add_argument("--initial-dir",type=Path,required=True)
    ap.add_argument("--candidates",type=Path,required=True)
    ap.add_argument("--execution-rule",type=Path,required=True)
    ap.add_argument("--output",type=Path,required=True)
    args=ap.parse_args()

    rule=json.loads(args.execution_rule.read_text())
    if rule.get("schema")!="ttf_relational_historical_occurrence_bounded_execution_v0.2":
        raise RuntimeError("unexpected bounded Study-C execution rule")
    if rule.get("status")!="FROZEN_FINAL_BOUNDED_OCCURRENCE_EXECUTION_BEFORE_RESULT":
        raise RuntimeError("bounded Study-C execution is not frozen")
    if any(bool(v) for v in rule["response_firewall"].values()):
        raise RuntimeError("bounded Study-C execution firewall is open")

    rows=list(csv.DictReader(args.candidates.open(encoding="utf-8")))
    names=[str(row["species"]).strip() for row in rows]
    if len(names)!=1000 or len(set(names))!=1000:
        raise RuntimeError("expected exact frozen Study-C 1000-species candidate set")

    ledgers=sorted(args.initial_dir.glob("ledger-bounded-*.json"))
    if len(ledgers)!=250:
        raise RuntimeError(f"expected 250 bounded initial ledgers, found {len(ledgers)}")

    seen={}
    batches=set()
    for path in ledgers:
        p=json.loads(path.read_text())
        if p.get("schema")!="ttf_relational_historical_occurrence_bounded_batch_v0.2":
            raise RuntimeError(f"unexpected bounded initial schema: {path}")
        if p.get("scientific_query_change") is not False:
            raise RuntimeError("bounded initial run changed scientific query")
        if any(bool(v) for v in p["response_firewall"].values()):
            raise RuntimeError("bounded initial response firewall is open")
        batch=int(p["batch_index"])
        if batch in batches:
            raise RuntimeError(f"duplicate bounded initial batch {batch}")
        batches.add(batch)
        expected=[name for i,name in enumerate(names) if i%250==batch]
        if list(map(str,p["species_requested"]))!=expected:
            raise RuntimeError(f"bounded initial batch assignment drift: {batch}")
        for row in p["species"]:
            name=str(row["species"])
            if name in seen:
                raise RuntimeError(f"duplicate bounded initial species: {name}")
            seen[name]=dict(row)

    if batches!=set(range(250)) or set(seen)!=set(names):
        raise RuntimeError("bounded initial execution does not cover exact frozen 1000 species")

    request_errors=[name for name in names if seen[name]["status"]=="REQUEST_ERROR"]
    batch_size=int(rule["retry_contract"]["batch_size"])
    if batch_size!=4:
        raise RuntimeError("bounded retry batch-size drift")
    retry_batches=int(math.ceil(len(request_errors)/batch_size)) if request_errors else 0
    if retry_batches>250:
        raise RuntimeError("bounded retry matrix exceeds frozen maximum")

    payload={
        "schema":"ttf_relational_historical_occurrence_retry_plan_v0.2",
        "status":"RETRY_EXACT_INITIAL_REQUEST_ERRORS" if request_errors else "NO_RETRY_NEEDED",
        "candidate_species":1000,
        "initial_batches":250,
        "initial_request_error_count":len(request_errors),
        "request_error_species":request_errors,
        "retry_batch_size":batch_size,
        "retry_batch_count":retry_batches,
        "maximum_retry_rounds":1,
        "additional_retry_authorized":False,
        "scientific_query_change":False,
        "relation_result_seen":False,
        "genetic_response_used":False,
        "response_firewall":{
            "Study_C_sequence_identity_opened":False,
            "Study_C_pairwise_genetic_distances_opened":False,
            "Study_C_T_st_computed":False,
            "Study_C_beta_hist_computed":False,
        },
    }
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(payload,indent=2,sort_keys=True)+"\n")
    print(json.dumps({
        "status":payload["status"],
        "initial_request_error_count":len(request_errors),
        "retry_batch_count":retry_batches,
    },sort_keys=True))
    return 0


if __name__=="__main__":
    raise SystemExit(main())
