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
    if census.get("schema")!="ttf_genetic_palearctic_lgm_nonlepidoptera_census_v0.2":
        raise RuntimeError("unexpected v0.2 census")
    rule=json.loads(args.external_rule.read_text())
    if rule.get("schema")!="ttf_genetic_palearctic_lgm_external_data_rule_v0.2":
        raise RuntimeError("unexpected v0.2 external-data rule")
    if any(bool(v) for v in census["response_firewall"].values()):
        raise RuntimeError("v0.2 census genetic response firewall is open")

    formal={str(row["species"]) for row in census["retained_species"]}
    if len(formal)!=23:
        raise RuntimeError("formal v0.2 census species drift")

    ledgers={}
    rows=[]
    for path in sorted(args.input_dir.glob("*.json")):
        payload=json.loads(path.read_text())
        if payload.get("schema")!="ttf_genetic_palearctic_lgm_occurrence_shard_v0.2":
            continue
        if any(bool(v) for v in payload["response_firewall"].values()):
            raise RuntimeError("v0.2 occurrence shard response firewall is open")
        for row in payload["species"]:
            sp=str(row["species"])
            if sp in ledgers:
                raise RuntimeError(f"duplicate species ledger: {sp}")
            ledgers[sp]=row

    for path in sorted(args.input_dir.glob("*.csv")):
        with path.open(newline="",encoding="utf-8") as fh:
            rows.extend(dict(row) for row in csv.DictReader(fh))

    if set(ledgers)!=formal:
        raise RuntimeError(
            f"GBIF ledger coverage drift: missing={sorted(formal-set(ledgers))}, "
            f"extra={sorted(set(ledgers)-formal)}"
        )

    request_error=sorted(sp for sp,x in ledgers.items() if x["status"]=="REQUEST_ERROR")
    if request_error:
        status="INCOMPLETE_TECHNICAL_V02_GBIF_REQUEST_ERRORS"
        admissible=[]
    else:
        admissible=sorted(
            sp for sp,x in ledgers.items()
            if x["status"]=="PASS_OCCURRENCE_GEOMETRY"
        )
        required=int(rule["missingness"]["after_drop_require"]["minimum_total_species"])
        status=(
            "PASS_TO_V02_CURRENT_LGM_CLIMATE_EXTRACTION"
            if len(admissible)>=required
            else str(rule["missingness"]["failure"])
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
        writer.writeheader()
        writer.writerows(retained_rows)

    counts={
        sp:int(ledgers[sp].get("retained_palearctic_after_exact_dedup_thinning_cap",0))
        for sp in sorted(formal)
    }
    n=len(admissible)
    payload={
        "schema":"ttf_genetic_palearctic_lgm_occurrence_result_v0.2",
        "status":status,
        "formal_species":23,
        "admissible_species":n,
        "source_clusters":n,
        "target_clusters":n,
        "directed_nonself_dyads":n*(n-1),
        "failed_occurrence_species":sorted(
            sp for sp,x in ledgers.items()
            if x["status"]=="FAIL_OCCURRENCE_GEOMETRY"
        ),
        "rejected_taxon_match_species":sorted(
            sp for sp,x in ledgers.items()
            if x["status"]=="REJECTED_GBIF_TAXON_MATCH"
        ),
        "request_error_species":request_error,
        "admissible_species_list":admissible,
        "retained_counts":counts,
        "retained_occurrence_rows":len(retained_rows),
        "external_rule_sha256":sha256_path(args.external_rule),
        "census_sha256":sha256_path(args.census),
        "output_csv_sha256":sha256_path(args.output_csv),
        "response_firewall":{
            "species_level_genetic_scores_used":False,
            "pairwise_v02_T_st_computed":False,
            "beta_LGM_computed":False,
        },
        "next_step":(
            "Sample current and 21-ka CHELSA climate and build the v0.2 all-role LGM relation."
            if status=="PASS_TO_V02_CURRENT_LGM_CLIMATE_EXTRACTION"
            else (
                "Retry REQUEST_ERROR species only under the frozen v0.2 acquisition contract."
                if request_error else
                "STOP. Do not compute LGM relation or subgroup genetic response."
            )
        ),
    }
    args.output_summary.write_text(json.dumps(payload,indent=2,sort_keys=True)+"\n")
    print(json.dumps({
        "status":status,
        "admissible_species":n,
        "source_clusters":n,
        "target_clusters":n,
        "directed_nonself_dyads":n*(n-1),
        "request_errors":len(request_error),
        "failed_occurrence_species":len(payload["failed_occurrence_species"]),
    },sort_keys=True))
    return 0


if __name__=="__main__":
    raise SystemExit(main())
