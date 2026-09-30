#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

try:
    from scripts.acquire_palearctic_lgm_occurrences import (
        fetch_species_palearctic,
        load_palearctic_geometry,
    )
except ModuleNotFoundError:
    from acquire_palearctic_lgm_occurrences import (
        fetch_species_palearctic,
        load_palearctic_geometry,
    )

from ttf.relational_external import shard_items


def main() -> int:
    ap=argparse.ArgumentParser()
    ap.add_argument("--census",type=Path,required=True)
    ap.add_argument("--realm-geojson",type=Path,required=True)
    ap.add_argument("--output-csv",type=Path,required=True)
    ap.add_argument("--output-ledger",type=Path,required=True)
    ap.add_argument("--shard-index",type=int,default=0)
    ap.add_argument("--shards",type=int,default=1)
    args=ap.parse_args()

    census=json.loads(args.census.read_text())
    if census.get("schema")!="ttf_genetic_palearctic_lgm_nonlepidoptera_census_v0.2":
        raise RuntimeError("unexpected non-Lepidoptera Palearctic-LGM census")
    if census.get("status")!="PASS_TO_FRESH_V02_GBIF_CENSUS":
        raise RuntimeError("v0.2 census did not authorize external occurrence acquisition")
    if any(bool(v) for v in census["response_firewall"].values()):
        raise RuntimeError("v0.2 genetic response firewall is open")

    retained=list(census["retained_species"])
    names=[str(row["species"]) for row in retained]
    if len(names)!=23 or len(set(names))!=23:
        raise RuntimeError("expected exact frozen 23-species non-Lepidoptera panel")

    shard=list(shard_items(names,shard_index=args.shard_index,shards=args.shards))
    realm=load_palearctic_geometry(args.realm_geojson)

    output=[]
    ledger=[]
    for species in shard:
        try:
            rows,meta=fetch_species_palearctic(species,realm)
        except Exception as exc:
            rows=[]
            meta={
                "species":species,
                "status":"REQUEST_ERROR",
                "error":f"{type(exc).__name__}: {exc}",
            }
        for row in rows:
            row["role"]="both"
        output.extend(rows)
        meta["role"]="both"
        ledger.append(meta)

    args.output_csv.parent.mkdir(parents=True,exist_ok=True)
    fields=[
        "species","role","source_key","latitude","longitude",
        "priority_rank","priority_sha256",
    ]
    with args.output_csv.open("w",newline="",encoding="utf-8") as fh:
        writer=csv.DictWriter(fh,fieldnames=fields)
        writer.writeheader()
        writer.writerows(output)

    payload={
        "schema":"ttf_genetic_palearctic_lgm_occurrence_shard_v0.2",
        "shard_index":int(args.shard_index),
        "shards":int(args.shards),
        "species_requested":len(shard),
        "retained_occurrence_rows":len(output),
        "species":ledger,
        "response_firewall":{
            "species_level_genetic_scores_used":False,
            "pairwise_v02_T_st_computed":False,
            "beta_LGM_computed":False,
        },
    }
    args.output_ledger.write_text(json.dumps(payload,indent=2,sort_keys=True)+"\n")
    print(json.dumps({
        "schema":payload["schema"],
        "shard_index":args.shard_index,
        "species_requested":len(shard),
        "retained_occurrence_rows":len(output),
        "request_errors":sum(x["status"]=="REQUEST_ERROR" for x in ledger),
    },sort_keys=True))
    return 0


if __name__=="__main__":
    raise SystemExit(main())
