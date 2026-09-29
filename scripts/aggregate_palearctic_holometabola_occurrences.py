#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from collections import Counter
from pathlib import Path


def sha256_path(path: Path) -> str:
    h=hashlib.sha256()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda:handle.read(1<<20),b""):
            h.update(chunk)
    return h.hexdigest()


def main() -> int:
    ap=argparse.ArgumentParser()
    ap.add_argument("--input-dir",type=Path,required=True)
    ap.add_argument("--support",type=Path,required=True)
    ap.add_argument("--rule",type=Path,required=True)
    ap.add_argument("--output-csv",type=Path,required=True)
    ap.add_argument("--output-json",type=Path,required=True)
    args=ap.parse_args()

    support=json.loads(args.support.read_text())
    rule=json.loads(args.rule.read_text())
    if support.get("status")!="PASS_GBIF_SUPPORT_TO_LGM_SDM":
        raise RuntimeError("support gate not passed")
    if rule.get("status")!="FROZEN_BEFORE_FULL_OCCURRENCE_ASSET":
        raise RuntimeError("occurrence rule not frozen")

    expected=set(support["usable_source_species"]+support["usable_target_species"])
    ledger_rows=[]; occurrence_rows=[]
    ledger_files=sorted(args.input_dir.glob("ledger_*.json"))
    csv_files=sorted(args.input_dir.glob("occ_*.csv"))
    if len(ledger_files)!=6 or len(csv_files)!=6:
        raise RuntimeError("expected exactly six occurrence shard ledgers and CSVs")

    for path in ledger_files:
        payload=json.loads(path.read_text())
        if payload.get("schema")!="ttf_genetic_palearctic_holometabola_occurrence_shard_v0.1":
            raise RuntimeError(f"unexpected shard schema: {path}")
        ledger_rows.extend(payload["species"])
    for path in csv_files:
        with path.open(newline="",encoding="utf-8") as handle:
            occurrence_rows.extend(csv.DictReader(handle))

    names=[str(row["species"]) for row in ledger_rows]
    if len(names)!=len(set(names)) or set(names)!=expected:
        raise RuntimeError("full occurrence asset species set drift")
    request_errors=sorted(row["species"] for row in ledger_rows if row["status"]=="REQUEST_ERROR")
    failed=sorted(row["species"] for row in ledger_rows if row["status"]!="PASS_FULL_OCCURRENCE_ASSET")
    if request_errors:
        status="INCOMPLETE_TECHNICAL_FULL_OCCURRENCE_ASSET"
    elif failed:
        status="STOP_UNEXPECTED_FULL_OCCURRENCE_ASSET_FAILURE"
    else:
        status="PASS_FULL_OCCURRENCE_ASSET"

    occurrence_rows.sort(key=lambda r:(r["species"],int(r["priority_rank"]),int(r["source_key"])))
    args.output_csv.parent.mkdir(parents=True,exist_ok=True)
    fields=["species","frozen_side","source_key","latitude","longitude","priority_rank","priority_sha256"]
    with args.output_csv.open("w",newline="",encoding="utf-8") as handle:
        writer=csv.DictWriter(handle,fieldnames=fields); writer.writeheader(); writer.writerows(occurrence_rows)

    per_species=Counter(row["species"] for row in occurrence_rows)
    payload={
        "schema":"ttf_genetic_palearctic_holometabola_occurrence_asset_v0.1",
        "status":status,
        "species":len(expected),
        "rows":len(occurrence_rows),
        "minimum_rows_per_species":min(per_species.values()) if per_species else 0,
        "maximum_rows_per_species":max(per_species.values()) if per_species else 0,
        "occurrence_csv_sha256":sha256_path(args.output_csv),
        "support_receipt_sha256":sha256_path(args.support),
        "acquisition_rule_sha256":sha256_path(args.rule),
        "request_error_species":request_errors,
        "failed_species":failed,
        "per_species_retained_counts":dict(sorted(per_species.items())),
        "ledger":sorted(ledger_rows,key=lambda r:r["species"]),
        "firewall":{
            "species_level_phase4_scores_read":False,
            "pairwise_genetic_concordance_opened":False,
            "lgm_relation_computed":False
        }
    }
    args.output_json.parent.mkdir(parents=True,exist_ok=True)
    args.output_json.write_text(json.dumps(payload,indent=2,sort_keys=True)+"\n")
    print(json.dumps({
        "status":status,"species":len(expected),"rows":len(occurrence_rows),
        "minimum_rows_per_species":payload["minimum_rows_per_species"],
        "maximum_rows_per_species":payload["maximum_rows_per_species"],
        "request_errors":len(request_errors),"failed_species":failed,
        "occurrence_csv_sha256":payload["occurrence_csv_sha256"],
    },sort_keys=True))
    if status!="PASS_FULL_OCCURRENCE_ASSET":
        return 2
    return 0


if __name__=="__main__":
    raise SystemExit(main())
