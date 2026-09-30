#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path


def sha256_path(path: Path) -> str:
    h=hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1<<20),b""):
            h.update(chunk)
    return h.hexdigest()


def main() -> int:
    ap=argparse.ArgumentParser()
    ap.add_argument("--input-dir",type=Path,required=True)
    ap.add_argument("--census",type=Path,required=True)
    ap.add_argument("--external-rule",type=Path,required=True)
    ap.add_argument("--output-csv",type=Path,required=True)
    ap.add_argument("--output-summary",type=Path,required=True)
    args=ap.parse_args()

    census=json.loads(args.census.read_text())
    rule=json.loads(args.external_rule.read_text())
    if rule.get("schema")!="ttf_genetic_palearctic_lgm_external_data_rule_v0.1":
        raise RuntimeError("unexpected external-data rule")
    formal={str(row["species"]):str(row["role"]) for row in census["retained_species"]}
    if len(formal)!=44:
        raise RuntimeError("formal census species drift")

    ledgers={}
    rows=[]
    for path in sorted(args.input_dir.glob("*.json")):
        payload=json.loads(path.read_text())
        if payload.get("schema")!="ttf_genetic_palearctic_lgm_occurrence_shard_v0.1":
            continue
        if any(bool(v) for v in payload["response_firewall"].values()):
            raise RuntimeError("shard response firewall is open")
        for row in payload["species"]:
            sp=str(row["species"])
            if sp in ledgers:
                raise RuntimeError(f"duplicate species ledger: {sp}")
            ledgers[sp]=row

    for path in sorted(args.input_dir.glob("*.csv")):
        with path.open(newline="",encoding="utf-8") as fh:
            for row in csv.DictReader(fh):
                rows.append(dict(row))

    if set(ledgers)!=set(formal):
        raise RuntimeError(
            f"GBIF ledger coverage drift: missing={sorted(set(formal)-set(ledgers))}, "
            f"extra={sorted(set(ledgers)-set(formal))}"
        )

    request_error=sorted(sp for sp,x in ledgers.items() if x["status"]=="REQUEST_ERROR")
    if request_error:
        status="INCOMPLETE_TECHNICAL_GBIF_REQUEST_ERRORS"
        admissible=[]
    else:
        admissible=sorted(
            sp for sp,x in ledgers.items()
            if x["status"]=="PASS_OCCURRENCE_GEOMETRY"
        )
        sources=[sp for sp in admissible if formal[sp]=="source"]
        targets=[sp for sp in admissible if formal[sp]=="target"]
        required=rule["missingness"]["after_drop_require"]
        passed=(
            len(admissible)>=int(required["minimum_total_species"])
            and len(sources)>=int(required["minimum_source_clusters"])
            and len(targets)>=int(required["minimum_target_clusters"])
        )
        status=(
            "PASS_TO_CURRENT_LGM_CLIMATE_EXTRACTION"
            if passed else str(rule["missingness"]["failure"])
        )

    admissible_set=set(admissible)
    retained_rows=[row for row in rows if row["species"] in admissible_set]
    retained_rows.sort(
        key=lambda x:(x["species"],int(x["priority_rank"]),int(x["source_key"]))
    )
    args.output_csv.parent.mkdir(parents=True,exist_ok=True)
    fields=[
        "species","role","source_key","latitude","longitude",
        "priority_rank","priority_sha256",
    ]
    with args.output_csv.open("w",newline="",encoding="utf-8") as fh:
        writer=csv.DictWriter(fh,fieldnames=fields)
        writer.writeheader(); writer.writerows(retained_rows)

    sources=[sp for sp in admissible if formal[sp]=="source"]
    targets=[sp for sp in admissible if formal[sp]=="target"]
    counts={sp:int(ledgers[sp].get(
        "retained_palearctic_after_exact_dedup_thinning_cap",0
    )) for sp in formal}
    payload={
        "schema":"ttf_genetic_palearctic_lgm_occurrence_result_v0.1",
        "status":status,
        "formal_species":44,
        "admissible_species":len(admissible),
        "source_clusters":len(sources),
        "target_clusters":len(targets),
        "failed_occurrence_species":sorted(
            sp for sp,x in ledgers.items()
            if x["status"]=="FAIL_OCCURRENCE_GEOMETRY"
        ),
        "request_error_species":request_error,
        "admissible_species_list":admissible,
        "source_species":sources,
        "target_species":targets,
        "retained_counts":counts,
        "retained_occurrence_rows":len(retained_rows),
        "external_rule_sha256":sha256_path(args.external_rule),
        "census_sha256":sha256_path(args.census),
        "output_csv_sha256":sha256_path(args.output_csv),
        "response_firewall":{
            "species_level_genetic_scores_used":False,
            "pairwise_subpanel_T_st_computed":False,
            "beta_LGM_computed":False,
        },
        "next_step":(
            "Sample current and 21-ka CHELSA climate on frozen GBIF occurrences and Palearctic grid."
            if status=="PASS_TO_CURRENT_LGM_CLIMATE_EXTRACTION"
            else (
                "Retry REQUEST_ERROR species only under the same frozen acquisition contract."
                if request_error else
                "STOP. Do not compute subgroup genetic response."
            )
        ),
    }
    args.output_summary.write_text(json.dumps(payload,indent=2,sort_keys=True)+"\n")
    print(json.dumps({
        "status":status,
        "admissible_species":len(admissible),
        "source_clusters":len(sources),
        "target_clusters":len(targets),
        "request_errors":len(request_error),
        "failed_occurrence_species":len(payload["failed_occurrence_species"]),
    },sort_keys=True))


if __name__=="__main__":
    raise SystemExit(main())
