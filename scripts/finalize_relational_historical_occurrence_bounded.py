#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import json
from collections import Counter
from pathlib import Path

FIELDS = ["species", "source_key", "latitude", "longitude", "priority_rank", "priority_sha256"]


def load_initial(input_dir: Path, names: list[str]):
    ledgers = sorted(input_dir.glob("ledger-bounded-*.json"))
    csvs = {p.stem.removeprefix("occurrences-bounded-"): p for p in input_dir.glob("occurrences-bounded-*.csv")}
    if len(ledgers) != 250:
        raise RuntimeError(f"expected 250 initial ledgers, found {len(ledgers)}")
    species: dict[str, dict] = {}
    rows: dict[str, list[dict[str, str]]] = {}
    batches: set[int] = set()
    for path in ledgers:
        p=json.loads(path.read_text())
        if p.get("schema")!="ttf_relational_historical_occurrence_bounded_batch_v0.2":
            raise RuntimeError(f"unexpected initial schema: {path}")
        if p.get("scientific_query_change") is not False or any(bool(v) for v in p["response_firewall"].values()):
            raise RuntimeError("initial bounded transport contract violation")
        batch=int(p["batch_index"])
        if batch in batches:
            raise RuntimeError(f"duplicate initial batch {batch}")
        batches.add(batch)
        expected=[name for i,name in enumerate(names) if i % 250 == batch]
        if list(map(str,p["species_requested"])) != expected:
            raise RuntimeError(f"initial batch assignment drift {batch}")
        csv_path=csvs.get(str(batch))
        if csv_path is None:
            raise RuntimeError(f"missing initial occurrence CSV {batch}")
        by_name: dict[str,list[dict[str,str]]] = {}
        for row in csv.DictReader(csv_path.open(encoding="utf-8")):
            by_name.setdefault(str(row["species"]),[]).append(dict(row))
        for row in p["species"]:
            name=str(row["species"])
            if name in species:
                raise RuntimeError(f"duplicate initial species {name}")
            species[name]=dict(row)
            rows[name]=by_name.get(name,[])
    if batches != set(range(250)) or set(species) != set(names):
        raise RuntimeError("initial bounded transport coverage drift")
    return species,rows


def load_retry(input_dir: Path):
    ledgers=sorted(input_dir.glob("ledger-retry-*.json"))
    csvs={p.stem.removeprefix("occurrences-retry-"):p for p in input_dir.glob("occurrences-retry-*.csv")}
    species: dict[str,dict]={}
    rows: dict[str,list[dict[str,str]]]={}
    for path in ledgers:
        p=json.loads(path.read_text())
        if p.get("schema")!="ttf_relational_historical_occurrence_bounded_retry_batch_v0.2":
            raise RuntimeError(f"unexpected retry schema: {path}")
        if int(p.get("retry_round",-1))!=1:
            raise RuntimeError("unexpected Study-C retry round")
        if p.get("scientific_query_change") is not False or any(bool(v) for v in p["response_firewall"].values()):
            raise RuntimeError("retry bounded transport contract violation")
        batch=int(p["batch_index"])
        csv_path=csvs.get(str(batch))
        if csv_path is None:
            raise RuntimeError(f"missing retry occurrence CSV {batch}")
        by_name: dict[str,list[dict[str,str]]]={}
        for row in csv.DictReader(csv_path.open(encoding="utf-8")):
            by_name.setdefault(str(row["species"]),[]).append(dict(row))
        for row in p["species"]:
            name=str(row["species"])
            if name in species:
                raise RuntimeError(f"duplicate retry species {name}")
            species[name]=dict(row)
            rows[name]=by_name.get(name,[])
    return species,rows


def main() -> int:
    ap=argparse.ArgumentParser()
    ap.add_argument("--initial-dir",type=Path,required=True)
    ap.add_argument("--retry-dir",type=Path,required=True)
    ap.add_argument("--retry-plan",type=Path,required=True)
    ap.add_argument("--candidates",type=Path,required=True)
    ap.add_argument("--output-csv",type=Path,required=True)
    ap.add_argument("--output-ledger",type=Path,required=True)
    args=ap.parse_args()

    candidate_rows=list(csv.DictReader(args.candidates.open(encoding="utf-8")))
    names=[str(row["species"]).strip() for row in candidate_rows]
    if len(names)!=1000 or len(set(names))!=1000:
        raise RuntimeError("expected exact frozen Study-C 1000-species candidate set")

    plan=json.loads(args.retry_plan.read_text())
    if plan.get("schema")!="ttf_relational_historical_occurrence_retry_plan_v0.2":
        raise RuntimeError("unexpected Study-C retry plan")
    if int(plan.get("candidate_species",-1))!=1000 or int(plan.get("maximum_retry_rounds",-1))!=1:
        raise RuntimeError("Study-C retry-plan finite execution drift")
    if plan.get("scientific_query_change") is not False or any(bool(v) for v in plan["response_firewall"].values()):
        raise RuntimeError("Study-C retry-plan firewall drift")

    species,rows=load_initial(args.initial_dir,names)
    initial_errors=list(map(str,plan["request_error_species"]))
    actual_initial_errors=[name for name in names if species[name]["status"]=="REQUEST_ERROR"]
    if actual_initial_errors!=initial_errors:
        raise RuntimeError("Study-C initial REQUEST_ERROR set does not match frozen retry plan")

    retry_species,retry_rows=load_retry(args.retry_dir)
    if initial_errors:
        if set(retry_species)!=set(initial_errors):
            raise RuntimeError("Study-C retry species set differs from exact frozen initial REQUEST_ERROR set")
        for name in initial_errors:
            if species[name]["status"]!="REQUEST_ERROR":
                raise RuntimeError(f"attempt to retry non-REQUEST_ERROR species {name}")
            species[name]=retry_species[name]
            rows[name]=retry_rows.get(name,[])
    elif retry_species:
        raise RuntimeError("retry artifacts exist despite NO_RETRY_NEEDED")

    unresolved=[name for name in names if species[name]["status"]=="REQUEST_ERROR"]
    counts=Counter(str(species[name]["status"]) for name in names)
    passed=int(counts.get("PASS_OCCURRENCE_GEOMETRY",0))
    status=(
        "NOT_EVALUABLE_HISTORICAL_TECHNICAL_TRANSPORT"
        if unresolved else
        ("PASS_TO_HISTORICAL_ASSET_EXTRACTION" if passed>=500 else "NOT_EVALUABLE_HISTORICAL_OCCURRENCE_GEOMETRY")
    )

    all_rows=[]
    for name in names:
        all_rows.extend(rows.get(name,[]))
    all_rows.sort(key=lambda row:(row["species"],int(row["priority_rank"]),int(row["source_key"])))

    args.output_csv.parent.mkdir(parents=True,exist_ok=True)
    with args.output_csv.open("w",newline="",encoding="utf-8") as h:
        w=csv.DictWriter(h,fieldnames=FIELDS);w.writeheader();w.writerows(all_rows)

    payload={
        "schema":"ttf_relational_historical_occurrence_acquisition_v0.2",
        "status":status,
        "species":1000,
        "species_passing_ge_30":passed,
        "status_counts":dict(sorted(counts.items())),
        "retained_rows":len(all_rows),
        "initial_request_error_count":len(initial_errors),
        "final_request_error_count":len(unresolved),
        "unresolved_request_error_species":unresolved,
        "retry_rounds_completed":1 if initial_errors else 0,
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
    args.output_ledger.write_text(json.dumps(payload,indent=2,sort_keys=True)+"\n")
    print(json.dumps({
        "status":status,
        "species_passing_ge_30":passed,
        "initial_request_errors":len(initial_errors),
        "final_request_errors":len(unresolved),
        "status_counts":payload["status_counts"],
    },sort_keys=True))
    return 0


if __name__=="__main__":
    raise SystemExit(main())
