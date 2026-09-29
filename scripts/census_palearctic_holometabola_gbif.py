#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
from pathlib import Path

from scripts.census_relational_gbif_occurrence import fetch_species
from ttf.relational_external import shard_items


def main() -> int:
    ap=argparse.ArgumentParser()
    ap.add_argument("--panel",type=Path,required=True)
    ap.add_argument("--output",type=Path,required=True)
    ap.add_argument("--shard-index",type=int,required=True)
    ap.add_argument("--shards",type=int,required=True)
    args=ap.parse_args()

    panel=json.loads(args.panel.read_text())
    if panel.get("schema")!="ttf_genetic_palearctic_holometabola_panel_v0.1":
        raise RuntimeError("unexpected Palearctic panel schema")
    if panel.get("status")!="PASS_PANEL_SIZE_TO_LGM_SDM":
        raise RuntimeError("Palearctic panel did not pass the frozen size gate")
    if panel.get("genetic_pair_response_opened") is not False:
        raise RuntimeError("pairwise genetic response is already open")

    side={}
    for name in panel["source_species"]:
        side[name]="source"
    for name in panel["target_species"]:
        if name in side:
            raise RuntimeError("source/target overlap")
        side[name]="target"
    names=sorted(side)
    if len(names)!=44:
        raise RuntimeError("frozen Palearctic panel count drift")

    shard=list(shard_items(names,shard_index=args.shard_index,shards=args.shards))
    results=[]
    for i,name in enumerate(shard,start=1):
        try:
            row=fetch_species(
                name,
                year_range="1990,2026",
                thin_km=10.0,
                required_thinned=50,
                fetch_cap=3000,
            )
        except Exception as exc:
            row={
                "species":name,
                "status":"request_error",
                "error":f"{type(exc).__name__}: {exc}",
            }
        row["frozen_side"]=side[name]
        results.append(row)
        print(json.dumps({
            "shard":args.shard_index,
            "completed":i,
            "total":len(shard),
            "species":name,
            "status":row["status"],
            "thinned_records":row.get("thinned_records"),
        },sort_keys=True))

    payload={
        "schema":"ttf_genetic_palearctic_holometabola_gbif_census_shard_v0.1",
        "panel_status":panel["status"],
        "shard_index":args.shard_index,
        "shards":args.shards,
        "year_range":"1990,2026",
        "thin_km":10.0,
        "required_thinned":50,
        "fetch_cap":3000,
        "results":results,
        "firewall":{
            "species_level_phase4_scores_read":False,
            "pairwise_genetic_concordance_opened":False,
            "lgm_sdm_fitted":False
        }
    }
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(payload,indent=2,sort_keys=True)+"\n")
    return 0


if __name__=="__main__":
    raise SystemExit(main())
