#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

from shapely.geometry import Point, shape
from shapely.prepared import prep

try:
    from scripts.acquire_relational_environment_occurrences import get_json
except ModuleNotFoundError:
    from acquire_relational_environment_occurrences import get_json

from ttf.relational_environment import (
    deterministic_page_offsets,
    filter_and_thin_occurrences,
)
from ttf.relational_external import accepted_gbif_species_match, shard_items


def load_palearctic_geometry(path: Path):
    payload=json.loads(path.read_text())
    if payload.get("type")=="FeatureCollection":
        features=payload.get("features") or []
        if len(features)!=1:
            raise RuntimeError("Palearctic union GeoJSON must contain exactly one feature")
        geometry=shape(features[0]["geometry"])
    elif payload.get("type")=="Feature":
        geometry=shape(payload["geometry"])
    else:
        geometry=shape(payload)
    if geometry.is_empty:
        raise RuntimeError("Palearctic union geometry is empty")
    return prep(geometry)


def in_palearctic(prepared, lat: float, lon: float) -> bool:
    point=Point(float(lon),float(lat))
    return bool(prepared.intersects(point))


def fetch_species_palearctic(
    species: str,
    prepared_realm,
    *,
    minimum_retained: int=30,
    maximum_retained: int=200,
) -> tuple[list[dict[str,object]],dict[str,object]]:
    match=get_json("species/match",{"name":species,"strict":"true"})
    accepted,match_rule=accepted_gbif_species_match(species,match)
    usage_key=int(match.get("usageKey") or 0)
    if not accepted or usage_key<=0:
        return [],{
            "species":species,
            "status":"REJECTED_GBIF_TAXON_MATCH",
            "gbif_usage_key":usage_key,
            "gbif_match_type":match.get("matchType"),
            "gbif_canonical_name":match.get("canonicalName"),
            "gbif_rank":match.get("rank"),
            "match_rule":match_rule,
        }

    count_payload=get_json(
        "occurrence/search",
        {
            "taxonKey":usage_key,
            "hasCoordinate":"true",
            "hasGeospatialIssue":"false",
            "occurrenceStatus":"PRESENT",
            "year":"2010,2026",
            "limit":1,
        },
    )
    total=int(count_payload.get("count") or 0)
    offsets=deterministic_page_offsets(total,page_size=300,maximum_pages=20)
    global_records=[]
    palearctic_records=[]
    for offset in offsets:
        payload=get_json(
            "occurrence/search",
            {
                "taxonKey":usage_key,
                "hasCoordinate":"true",
                "hasGeospatialIssue":"false",
                "occurrenceStatus":"PRESENT",
                "year":"2010,2026",
                "limit":300,
                "offset":int(offset),
            },
        )
        for item in payload.get("results") or []:
            lat=item.get("decimalLatitude")
            lon=item.get("decimalLongitude")
            key=int(item.get("key") or 0)
            if lat is None or lon is None or key<=0:
                continue
            try:
                lat=float(lat); lon=float(lon)
            except (TypeError,ValueError):
                continue
            row={"source_key":key,"latitude":lat,"longitude":lon}
            global_records.append(row)
            if in_palearctic(prepared_realm,lat,lon):
                palearctic_records.append(row)

    retained=filter_and_thin_occurrences(
        species,
        palearctic_records,
        minimum_distance_km=10.0,
        maximum_retained=maximum_retained,
    )
    status=(
        "PASS_OCCURRENCE_GEOMETRY"
        if len(retained)>=minimum_retained
        else "FAIL_OCCURRENCE_GEOMETRY"
    )
    return retained,{
        "species":species,
        "status":status,
        "gbif_usage_key":usage_key,
        "gbif_total_coordinate_count_2010_2026":total,
        "queried_page_offsets":list(offsets),
        "queried_rows_with_coordinates":len(global_records),
        "queried_rows_in_palearctic":len(palearctic_records),
        "retained_palearctic_after_exact_dedup_thinning_cap":len(retained),
        "gbif_match_type":match.get("matchType"),
        "gbif_canonical_name":match.get("canonicalName"),
        "gbif_rank":match.get("rank"),
        "match_rule":match_rule,
    }


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
    if census.get("schema")!="ttf_genetic_palearctic_lgm_subpanel_census_v0.1":
        raise RuntimeError("unexpected Palearctic-LGM census")
    if census.get("status")!="PASS_TO_RESPONSE_BLIND_LGM_EXTERNAL_DATA_CENSUS":
        raise RuntimeError("formal Palearctic census did not authorize external data")
    if any(bool(v) for v in census["response_firewall"].values()):
        raise RuntimeError("subpanel response firewall is open")

    retained=census["retained_species"]
    names=[str(row["species"]) for row in retained]
    roles={str(row["species"]):str(row["role"]) for row in retained}
    if len(names)!=44 or len(set(names))!=44:
        raise RuntimeError("expected exact frozen 44-species formal panel")

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
            row["role"]=roles[species]
        output.extend(rows)
        meta["role"]=roles[species]
        ledger.append(meta)

    args.output_csv.parent.mkdir(parents=True,exist_ok=True)
    fields=[
        "species","role","source_key","latitude","longitude",
        "priority_rank","priority_sha256",
    ]
    with args.output_csv.open("w",newline="",encoding="utf-8") as fh:
        writer=csv.DictWriter(fh,fieldnames=fields)
        writer.writeheader(); writer.writerows(output)

    payload={
        "schema":"ttf_genetic_palearctic_lgm_occurrence_shard_v0.1",
        "shard_index":args.shard_index,
        "shards":args.shards,
        "species_requested":len(shard),
        "retained_occurrence_rows":len(output),
        "species":ledger,
        "response_firewall":{
            "species_level_genetic_scores_used":False,
            "pairwise_subpanel_T_st_computed":False,
            "beta_LGM_computed":False,
        },
    }
    args.output_ledger.write_text(json.dumps(payload,indent=2,sort_keys=True)+"\n")


if __name__=="__main__":
    raise SystemExit(main())
