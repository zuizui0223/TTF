#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
from collections import Counter
from pathlib import Path


def main() -> int:
    ap=argparse.ArgumentParser()
    ap.add_argument("--input-dir",type=Path,required=True)
    ap.add_argument("--contract",type=Path,required=True)
    ap.add_argument("--panel",type=Path,required=True)
    ap.add_argument("--output",type=Path,required=True)
    args=ap.parse_args()

    contract=json.loads(args.contract.read_text())
    panel=json.loads(args.panel.read_text())
    rows=[]
    shards=sorted(args.input_dir.glob("*.json"))
    if not shards:
        raise RuntimeError("no GBIF census shards found")
    for path in shards:
        payload=json.loads(path.read_text())
        if payload.get("schema")!="ttf_genetic_palearctic_holometabola_gbif_census_shard_v0.1":
            continue
        rows.extend(payload["results"])

    expected=set(panel["source_species"]+panel["target_species"])
    names=[str(r["species"]) for r in rows]
    if len(names)!=len(set(names)) or set(names)!=expected:
        raise RuntimeError("GBIF census species set drift")

    request_errors=[r for r in rows if r["status"]=="request_error"]
    usable=[
        r for r in rows
        if r["status"]=="usable_lower_bound" and int(r.get("thinned_records",0))>=50
    ]
    usable_source=sorted(r["species"] for r in usable if r["frozen_side"]=="source")
    usable_target=sorted(r["species"] for r in usable if r["frozen_side"]=="target")
    min_total=int(contract["lgm_sdm"]["minimum_species_with_usable_sdm"])
    min_source=int(contract["lgm_sdm"]["minimum_source_species_with_usable_sdm"])
    min_target=int(contract["lgm_sdm"]["minimum_target_species_with_usable_sdm"])

    if request_errors:
        status="INCOMPLETE_TECHNICAL_GBIF_CENSUS"
    elif (
        len(usable)>=min_total
        and len(usable_source)>=min_source
        and len(usable_target)>=min_target
    ):
        status="PASS_GBIF_SUPPORT_TO_LGM_SDM"
    else:
        status="STOP_NOT_EVALUABLE_LGM_SDM_SUPPORT"

    payload={
        "schema":"ttf_genetic_palearctic_holometabola_gbif_census_v0.1",
        "status":status,
        "panel_species":len(expected),
        "usable_species":len(usable),
        "usable_source_species":usable_source,
        "usable_target_species":usable_target,
        "usable_source_count":len(usable_source),
        "usable_target_count":len(usable_target),
        "usable_primary_dyads":len(usable_source)*len(usable_target),
        "status_counts":dict(sorted(Counter(str(r["status"]) for r in rows).items())),
        "minimum_total":min_total,
        "minimum_source":min_source,
        "minimum_target":min_target,
        "request_error_species":sorted(r["species"] for r in request_errors),
        "results":sorted(rows,key=lambda r:r["species"]),
        "firewall":{
            "species_level_phase4_scores_read":False,
            "pairwise_genetic_concordance_opened":False,
            "lgm_sdm_fitted":False
        },
        "interpretation":"GBIF occurrence support is a response-blind prerequisite. Technical request errors never count as ecological insufficiency."
    }
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(payload,indent=2,sort_keys=True)+"\n")
    print(json.dumps({
        "status":status,
        "usable_species":len(usable),
        "usable_source_count":len(usable_source),
        "usable_target_count":len(usable_target),
        "usable_primary_dyads":len(usable_source)*len(usable_target),
        "request_errors":len(request_errors),
        "status_counts":payload["status_counts"],
    },sort_keys=True))
    return 0


if __name__=="__main__":
    raise SystemExit(main())
