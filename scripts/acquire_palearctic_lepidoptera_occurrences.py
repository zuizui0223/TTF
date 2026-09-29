#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import json
import runpy
from pathlib import Path


ROOT=Path(__file__).resolve().parents[1]
SOURCE=ROOT/"scripts/acquire_relational_environment_occurrences.py"


def main() -> int:
    ap=argparse.ArgumentParser()
    ap.add_argument("--screen",type=Path,required=True)
    ap.add_argument("--output-csv",type=Path,required=True)
    ap.add_argument("--output-ledger",type=Path,required=True)
    ap.add_argument("--shard-index",type=int,required=True)
    ap.add_argument("--shards",type=int,required=True)
    args=ap.parse_args()

    screen=json.loads(args.screen.read_text())
    if screen.get("schema")!="ttf_palearctic_lepidoptera_response_blind_screen_v0.1":
        raise RuntimeError("unexpected screen receipt")
    if screen.get("status")!="PASS_TO_LGM_PREDICTOR_CONSTRUCTION":
        raise RuntimeError("screen did not authorize LGM predictor construction")
    names=[str(x["species"]) for x in screen["primary_species"]]
    if len(names)!=26 or len(set(names))!=26:
        raise RuntimeError("expected exact frozen 26-species primary panel")

    fetch_species=runpy.run_path(str(SOURCE))["fetch_species"]
    shard=[name for i,name in enumerate(names) if i % args.shards == args.shard_index]

    retained_rows=[]
    ledgers=[]
    for i,species in enumerate(shard,1):
        try:
            retained,ledger=fetch_species(species)
        except Exception as exc:
            retained,ledger=[],{
                "species":species,
                "status":"REQUEST_ERROR",
                "error":f"{type(exc).__name__}: {exc}",
            }
        retained_rows.extend(retained)
        ledgers.append(ledger)
        print(json.dumps({
            "shard":args.shard_index,
            "completed":i,
            "total":len(shard),
            "species":species,
            "status":ledger.get("status"),
            "retained":ledger.get("retained_after_exact_dedup_thinning_cap",0),
        },sort_keys=True))

    args.output_csv.parent.mkdir(parents=True,exist_ok=True)
    fields=["species","source_key","latitude","longitude","priority_rank","priority_sha256"]
    with args.output_csv.open("w",newline="",encoding="utf-8") as handle:
        writer=csv.DictWriter(handle,fieldnames=fields)
        writer.writeheader()
        writer.writerows(retained_rows)

    payload={
        "schema":"ttf_palearctic_lepidoptera_occurrence_shard_v0.1",
        "shard_index":args.shard_index,
        "shards":args.shards,
        "species_requested":len(shard),
        "species_names":shard,
        "retained_occurrence_rows":len(retained_rows),
        "species":ledgers,
        "response_firewall":{
            "subpanel_nucleotide_identity_read":False,
            "subpanel_pairwise_genetic_distance_read":False,
            "subpanel_transfer_response_constructed":False,
        },
    }
    args.output_ledger.write_text(json.dumps(payload,indent=2,sort_keys=True)+"\n")
    return 0


if __name__=="__main__":
    raise SystemExit(main())
