#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path

from shapely.geometry import Point, shape
from shapely.prepared import prep

from ttf.gbif_transport import get_json
from ttf.relational_environment import deterministic_page_offsets, filter_and_thin_occurrences
from ttf.relational_external import accepted_gbif_species_match, shard_items


def sha256_path(path: Path) -> str:
    h=hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda:f.read(1<<20),b""):
            h.update(chunk)
    return h.hexdigest()


def load_native_ranges(path: Path):
    payload=json.loads(path.read_text())
    if payload.get("type")!="FeatureCollection":
        raise RuntimeError("native range asset must be GeoJSON FeatureCollection")
    out={}
    for feat in payload.get("features") or []:
        props=feat.get("properties") or {}
        host_id=str(props.get("accepted_host_id") or "").strip()
        if not host_id:
            raise RuntimeError("native range feature lacks accepted_host_id")
        geom=shape(feat["geometry"])
        if geom.is_empty:
            raise RuntimeError(f"empty native range geometry: {host_id}")
        out[host_id]=prep(geom)
    return out


def panel_hosts(path: Path):
    rows=list(csv.DictReader(path.open(encoding="utf-8")))
    hosts={}
    for row in rows:
        ids=[x.strip() for x in row["accepted_host_ids"].split(";") if x.strip()]
        names=[x.strip() for x in row["accepted_host_names"].split(";") if x.strip()]
        if len(ids)!=len(names) or int(row["n_hosts"])!=len(ids):
            raise RuntimeError(f"host panel drift: {row['species']}")
        for host_id,name in zip(ids,names):
            old=hosts.setdefault(host_id,name)
            if old!=name:
                raise RuntimeError(f"host ID/name inconsistency: {host_id}")
    return rows,hosts


def fetch_host(name: str, host_id: str, prepared, rule: dict):
    gbif=rule["gbif"]
    match=get_json("species/match",{"name":name,"strict":"true"})
    accepted,match_rule=accepted_gbif_species_match(name,match)
    usage_key=int(match.get("usageKey") or 0)
    if not accepted or usage_key<=0:
        return [],{
            "accepted_host_id":host_id,"accepted_host_name":name,
            "status":"REJECTED_GBIF_TAXON_MATCH",
            "gbif_usage_key":usage_key,
            "gbif_match_type":match.get("matchType"),
            "gbif_canonical_name":match.get("canonicalName"),
            "gbif_rank":match.get("rank"),
            "match_rule":match_rule,
        }

    filters=gbif["occurrence_filters"]
    count=get_json("occurrence/search",{
        "taxonKey":usage_key,
        "hasCoordinate":str(bool(filters["hasCoordinate"])).lower(),
        "hasGeospatialIssue":str(bool(filters["hasGeospatialIssue"])).lower(),
        "occurrenceStatus":filters["occurrenceStatus"],
        "year":filters["year"],
        "limit":1,
    })
    total=int(count.get("count") or 0)
    page_size=int(gbif["page_size"])
    offsets=deterministic_page_offsets(
        total,page_size=page_size,maximum_pages=int(gbif["maximum_pages"])
    )
    global_rows=[]; native_rows=[]
    for offset in offsets:
        payload=get_json("occurrence/search",{
            "taxonKey":usage_key,
            "hasCoordinate":str(bool(filters["hasCoordinate"])).lower(),
            "hasGeospatialIssue":str(bool(filters["hasGeospatialIssue"])).lower(),
            "occurrenceStatus":filters["occurrenceStatus"],
            "year":filters["year"],
            "limit":page_size,
            "offset":int(offset),
        })
        for item in payload.get("results") or []:
            lat=item.get("decimalLatitude"); lon=item.get("decimalLongitude")
            key=int(item.get("key") or 0)
            if lat is None or lon is None or key<=0:
                continue
            try:
                lat=float(lat); lon=float(lon)
            except (TypeError,ValueError):
                continue
            row={"source_key":key,"latitude":lat,"longitude":lon}
            global_rows.append(row)
            if prepared.intersects(Point(lon,lat)):
                native_rows.append(row)

    retained=filter_and_thin_occurrences(
        name,native_rows,
        minimum_distance_km=float(gbif["minimum_distance_km"]),
        maximum_retained=int(gbif["maximum_retained_per_host"]),
    )
    minimum=int(gbif["minimum_retained_per_host"])
    status="PASS_HOST_OCCURRENCE" if len(retained)>=minimum else "FAIL_HOST_OCCURRENCE"
    for row in retained:
        row["accepted_host_id"]=host_id
        row["accepted_host_name"]=name
    return retained,{
        "accepted_host_id":host_id,"accepted_host_name":name,
        "status":status,
        "gbif_usage_key":usage_key,
        "gbif_total_coordinate_count":total,
        "queried_page_offsets":list(offsets),
        "queried_rows_with_coordinates":len(global_rows),
        "queried_rows_in_primary_native_range":len(native_rows),
        "retained_after_exact_dedup_10km_thinning_cap":len(retained),
        "gbif_match_type":match.get("matchType"),
        "gbif_canonical_name":match.get("canonicalName"),
        "gbif_rank":match.get("rank"),
        "match_rule":match_rule,
    }


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--panel",type=Path,required=True)
    ap.add_argument("--rule",type=Path,required=True)
    ap.add_argument("--native-ranges",type=Path,required=True)
    ap.add_argument("--output-csv",type=Path,required=True)
    ap.add_argument("--output-ledger",type=Path,required=True)
    ap.add_argument("--shard-index",type=int,required=True)
    ap.add_argument("--shards",type=int,required=True)
    args=ap.parse_args()

    rule=json.loads(args.rule.read_text())
    if rule.get("schema")!="ttf_genetic_historical_host_occurrence_rule_v0.1":
        raise RuntimeError("unexpected historical host occurrence rule")
    if rule.get("status")!="FROZEN_BEFORE_WCVP_PANEL_RESULT_OR_ANY_HOST_GBIF_RESULT":
        raise RuntimeError("host occurrence rule not prospectively frozen")
    panel_rows,hosts=panel_hosts(args.panel)
    ranges=load_native_ranges(args.native_ranges)
    if set(hosts)!=set(ranges):
        raise RuntimeError("native range asset host set differs from frozen panel host set")

    items=sorted(hosts.items(),key=lambda x:(x[1],x[0]))
    shard=list(shard_items(items,shard_index=args.shard_index,shards=args.shards))
    output=[]; ledger=[]
    for i,(host_id,name) in enumerate(shard,start=1):
        try:
            rows,meta=fetch_host(name,host_id,ranges[host_id],rule)
        except Exception as exc:
            rows=[]; meta={
                "accepted_host_id":host_id,"accepted_host_name":name,
                "status":"REQUEST_ERROR","error":f"{type(exc).__name__}: {exc}",
            }
        output.extend(rows); ledger.append(meta)
        if i%5==0:
            print(json.dumps({"shard":args.shard_index,"completed":i,"total":len(shard)}))

    args.output_csv.parent.mkdir(parents=True,exist_ok=True)
    fields=["accepted_host_id","accepted_host_name","source_key","latitude","longitude","priority_rank","priority_sha256"]
    with args.output_csv.open("w",newline="",encoding="utf-8") as f:
        w=csv.DictWriter(f,fieldnames=fields); w.writeheader(); w.writerows(output)
    payload={
        "schema":"ttf_genetic_historical_host_occurrence_shard_v0.1",
        "shard_index":args.shard_index,"shards":args.shards,
        "panel_sha256":sha256_path(args.panel),
        "rule_sha256":sha256_path(args.rule),
        "native_ranges_sha256":sha256_path(args.native_ranges),
        "hosts_requested":len(shard),"retained_occurrence_rows":len(output),
        "hosts":ledger,
        "response_firewall":{
            "fresh_sequence_identity_opened":False,
            "fresh_pairwise_genetic_distances_opened":False,
            "fresh_post_ibd_turnover_computed":False,
            "fresh_beta_host_LGM_computed":False,
        },
    }
    args.output_ledger.write_text(json.dumps(payload,indent=2,sort_keys=True)+"\n")
    return 0


if __name__=="__main__":
    raise SystemExit(main())
