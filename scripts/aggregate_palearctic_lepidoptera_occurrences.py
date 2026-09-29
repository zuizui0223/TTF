#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path


def main() -> int:
    ap=argparse.ArgumentParser()
    ap.add_argument("--screen",type=Path,required=True)
    ap.add_argument("--contract",type=Path,required=True)
    ap.add_argument("--input-dir",type=Path,required=True)
    ap.add_argument("--output-csv",type=Path,required=True)
    ap.add_argument("--output-json",type=Path,required=True)
    args=ap.parse_args()

    screen=json.loads(args.screen.read_text())
    contract=json.loads(args.contract.read_text())
    names=[str(x["species"]) for x in screen["primary_species"]]
    expected=set(names)

    ledgers=[]
    for path in sorted(args.input_dir.glob("ledger_*.json")):
        payload=json.loads(path.read_text())
        if payload.get("schema")!="ttf_palearctic_lepidoptera_occurrence_shard_v0.1":
            continue
        ledgers.extend(payload["species"])
    by_name={str(x["species"]):x for x in ledgers}
    if set(by_name)!=expected:
        raise RuntimeError(
            f"occurrence shard coverage drift: missing={sorted(expected-set(by_name))} "
            f"extra={sorted(set(by_name)-expected)}"
        )

    rows=[]
    for path in sorted(args.input_dir.glob("occurrences_*.csv")):
        with path.open(newline="",encoding="utf-8") as handle:
            rows.extend(csv.DictReader(handle))

    statuses={}
    for name in names:
        status=str(by_name[name].get("status"))
        statuses[status]=statuses.get(status,0)+1

    passing=[
        name for name in names
        if by_name[name].get("status")=="PASS_OCCURRENCE_GEOMETRY"
    ]
    request_errors=[
        name for name in names
        if by_name[name].get("status")=="REQUEST_ERROR"
    ]
    minimum=int(contract["ecological_occurrences"]["minimum_predictor_admissible_species"])
    if request_errors:
        status="INCOMPLETE_TECHNICAL_REQUEST_ERRORS"
    elif len(passing)<minimum:
        status="STOP_LGM_SDM_PANEL_TOO_SMALL"
    else:
        status="PASS_TO_LGM_CLIMATE_EXTRACTION"

    args.output_csv.parent.mkdir(parents=True,exist_ok=True)
    fields=["species","source_key","latitude","longitude","priority_rank","priority_sha256"]
    with args.output_csv.open("w",newline="",encoding="utf-8") as handle:
        writer=csv.DictWriter(handle,fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)

    payload={
        "schema":"ttf_palearctic_lepidoptera_occurrence_execution_v0.1",
        "status":status,
        "primary_species":len(names),
        "species_status_counts":dict(sorted(statuses.items())),
        "predictor_admissible_species":len(passing),
        "predictor_admissible_species_names":passing,
        "request_error_species":request_errors,
        "minimum_predictor_admissible_species":minimum,
        "retained_occurrence_rows":len(rows),
        "response_firewall":{
            "subpanel_nucleotide_identity_read":False,
            "subpanel_pairwise_genetic_distance_read":False,
            "subpanel_transfer_response_constructed":False,
        },
        "next_step":(
            "extract frozen current and 21-ka CHELSA climate and fit the response-blind SDMs"
            if status=="PASS_TO_LGM_CLIMATE_EXTRACTION"
            else "STOP or repair only technical REQUEST_ERROR species under a separately frozen transport rule"
        )
    }
    args.output_json.write_text(json.dumps(payload,indent=2,sort_keys=True)+"\n")
    print(json.dumps({
        "status":status,
        "passing":len(passing),
        "request_errors":len(request_errors),
        "statuses":statuses,
        "rows":len(rows),
    },sort_keys=True))
    return 0


if __name__=="__main__":
    raise SystemExit(main())
