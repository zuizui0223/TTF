#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import time
from pathlib import Path
from urllib.error import HTTPError
from urllib.parse import urlencode
from urllib.request import Request, urlopen

from ttf.relational_environment import deterministic_page_offsets, haversine_km
from ttf.relational_external import accepted_gbif_species_match, shard_items

GBIF="https://api.gbif.org/v1"
USER_AGENT="ttf-palearctic-lgm/0.1"


def get_json(path: str, params: dict[str, object], retries: int = 8) -> dict:
    url=f"{GBIF}/{path}?{urlencode(params)}"
    error=None
    for attempt in range(retries):
        try:
            req=Request(url,headers={"User-Agent":USER_AGENT,"Accept":"application/json"})
            with urlopen(req,timeout=90) as response:
                payload=json.loads(response.read().decode("utf-8"))
            time.sleep(0.05)
            return payload
        except HTTPError as exc:
            error=exc
            if exc.code not in {429,500,502,503,504}:
                break
            retry_after=exc.headers.get("Retry-After")
            try:
                wait=float(retry_after) if retry_after else float(2**attempt)
            except Exception:
                wait=float(2**attempt)
            time.sleep(min(120.0,max(1.0,wait)))
        except Exception as exc:
            error=exc
            time.sleep(min(60.0,float(2**attempt)))
    raise RuntimeError(f"GBIF request failed: {url}") from error


def priority(species: str, source_key: int) -> str:
    return hashlib.sha256(
        f"palearctic-lgm-occ-v0.1|{species}|{int(source_key)}".encode()
    ).hexdigest()


def thin(species: str, records: list[dict[str, object]]) -> list[dict[str, object]]:
    by_coord={}
    for r in records:
        try:
            key=int(r["source_key"]); lat=float(r["latitude"]); lon=float(r["longitude"])
        except Exception:
            continue
        if key<=0 or not math.isfinite(lat) or not math.isfinite(lon):
            continue
        if not (-90<=lat<=90 and -180<=lon<=180) or (lat==0 and lon==0):
            continue
        coord=(lat,lon)
        current=by_coord.get(coord)
        if current is None or key<int(current["source_key"]):
            by_coord[coord]={"species":species,"source_key":key,"latitude":lat,"longitude":lon}
    ordered=sorted(
        by_coord.values(),
        key=lambda r:(priority(species,int(r["source_key"])),int(r["source_key"]))
    )
    kept=[]; coords=[]
    for row in ordered:
        coord=(float(row["latitude"]),float(row["longitude"]))
        if all(haversine_km(coord,other)>=10.0 for other in coords):
            out=dict(row)
            out["priority_rank"]=len(kept)
            out["priority_sha256"]=priority(species,int(row["source_key"]))
            kept.append(out); coords.append(coord)
            if len(kept)>=1000:
                break
    return kept


def fetch_species(species: str) -> tuple[list[dict],dict]:
    match=get_json("species/match",{"name":species,"strict":"true"})
    accepted,match_rule=accepted_gbif_species_match(species,match)
    key=int(match.get("usageKey") or 0)
    if not accepted or key<=0:
        return [],{
            "species":species,"status":"REJECTED_GBIF_TAXON_MATCH",
            "gbif_usage_key":key,"match_rule":match_rule,
            "canonical_name":match.get("canonicalName"),"match_type":match.get("matchType")
        }
    params={
        "taxonKey":key,"hasCoordinate":"true","hasGeospatialIssue":"false",
        "occurrenceStatus":"PRESENT","year":"1990,2026","limit":1,
    }
    count=get_json("occurrence/search",params)
    total=int(count.get("count") or 0)
    offsets=deterministic_page_offsets(total,page_size=300,maximum_pages=20)
    raw=[]
    for offset in offsets:
        payload=get_json("occurrence/search",{
            "taxonKey":key,"hasCoordinate":"true","hasGeospatialIssue":"false",
            "occurrenceStatus":"PRESENT","year":"1990,2026",
            "limit":300,"offset":int(offset),
        })
        for item in payload.get("results") or []:
            lat=item.get("decimalLatitude"); lon=item.get("decimalLongitude")
            occ=int(item.get("key") or 0)
            if lat is None or lon is None or occ<=0:
                continue
            raw.append({"source_key":occ,"latitude":lat,"longitude":lon})
    retained=thin(species,raw)
    return retained,{
        "species":species,
        "status":"PASS_FULL_OCCURRENCE_ASSET" if len(retained)>=50 else "FAIL_FULL_OCCURRENCE_ASSET",
        "gbif_usage_key":key,
        "gbif_total_coordinate_count_1990_2026":total,
        "queried_page_offsets":list(offsets),
        "raw_coordinate_rows":len(raw),
        "retained_after_dedup_thinning_cap":len(retained),
        "retained_source_keys_sha256":hashlib.sha256(
            json.dumps(sorted(int(r["source_key"]) for r in retained),separators=(",",":")).encode()
        ).hexdigest(),
        "match_rule":match_rule,
        "canonical_name":match.get("canonicalName"),
        "match_type":match.get("matchType"),
    }


def main() -> int:
    ap=argparse.ArgumentParser()
    ap.add_argument("--support",type=Path,required=True)
    ap.add_argument("--rule",type=Path,required=True)
    ap.add_argument("--shard-index",type=int,required=True)
    ap.add_argument("--shards",type=int,required=True)
    ap.add_argument("--output-csv",type=Path,required=True)
    ap.add_argument("--output-json",type=Path,required=True)
    args=ap.parse_args()

    support=json.loads(args.support.read_text())
    rule=json.loads(args.rule.read_text())
    if support.get("status")!="PASS_GBIF_SUPPORT_TO_LGM_SDM":
        raise RuntimeError("GBIF support gate not passed")
    if rule.get("status")!="FROZEN_BEFORE_FULL_OCCURRENCE_ASSET":
        raise RuntimeError("occurrence acquisition rule not frozen")
    side={name:"source" for name in support["usable_source_species"]}
    side.update({name:"target" for name in support["usable_target_species"]})
    names=sorted(side)
    if len(names)!=42:
        raise RuntimeError("usable species count drift")
    shard=list(shard_items(names,shard_index=args.shard_index,shards=args.shards))
    rows=[]; ledgers=[]
    for i,name in enumerate(shard,start=1):
        try:
            retained,ledger=fetch_species(name)
        except Exception as exc:
            retained=[]; ledger={"species":name,"status":"REQUEST_ERROR","error":f"{type(exc).__name__}: {exc}"}
        for row in retained:
            row["frozen_side"]=side[name]
        ledger["frozen_side"]=side[name]
        rows.extend(retained); ledgers.append(ledger)
        print(json.dumps({
            "shard":args.shard_index,"completed":i,"total":len(shard),
            "species":name,"status":ledger["status"],
            "retained":ledger.get("retained_after_dedup_thinning_cap")
        },sort_keys=True))

    args.output_csv.parent.mkdir(parents=True,exist_ok=True)
    fields=["species","frozen_side","source_key","latitude","longitude","priority_rank","priority_sha256"]
    with args.output_csv.open("w",newline="",encoding="utf-8") as handle:
        writer=csv.DictWriter(handle,fieldnames=fields); writer.writeheader(); writer.writerows(rows)
    payload={
        "schema":"ttf_genetic_palearctic_holometabola_occurrence_shard_v0.1",
        "shard_index":args.shard_index,"shards":args.shards,
        "species_requested":len(shard),"retained_rows":len(rows),
        "species":ledgers,
        "firewall":{"species_level_phase4_scores_read":False,"pairwise_genetic_concordance_opened":False}
    }
    args.output_json.write_text(json.dumps(payload,indent=2,sort_keys=True)+"\n")
    return 0


if __name__=="__main__":
    raise SystemExit(main())
