#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path


EXPECTED_SURVIVOR_GEOMETRY_SHA256 = "6e9ec4a6c56ffac82e91976aee2580a0b9d03e2b932dd907429fd09802419622"
EXPECTED_AUTHORIZATION_SHA256 = "47ffa048052d7bea5e76ebc144dcd40ea7e86910033ae082acc58757ca92d4d9"


def sha256_path(path: Path) -> str:
    h=hashlib.sha256()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda:handle.read(1<<20),b""):
            h.update(chunk)
    return h.hexdigest()


def load_contract(path: Path) -> dict:
    payload=json.loads(path.read_text())
    if payload.get("schema")!="ttf_genetic_palearctic_holometabola_lgm_program_v0.1":
        raise RuntimeError("unexpected Palearctic LGM contract schema")
    if payload.get("status")!="FROZEN_BEFORE_PAIRWISE_GENETIC_RESPONSE_CONSTRUCTION":
        raise RuntimeError("Palearctic LGM contract is not frozen")
    return payload


def deterministic_split(species: list[str], namespace: str) -> tuple[list[str],list[str]]:
    keyed=sorted(
        (hashlib.sha256(f"{namespace}|{name}".encode()).hexdigest(),name)
        for name in species
    )
    n=len(keyed)//2
    return (
        sorted(name for _,name in keyed[:n]),
        sorted(name for _,name in keyed[n:]),
    )


def main() -> int:
    ap=argparse.ArgumentParser()
    ap.add_argument("--survivor-geometry-csv",type=Path,required=True)
    ap.add_argument("--authorization-json",type=Path,required=True)
    ap.add_argument("--ecoregions",type=Path,required=True)
    ap.add_argument("--realm-source-receipt",type=Path,required=True)
    ap.add_argument("--realm-subset-receipt",type=Path,required=True)
    ap.add_argument("--contract",type=Path,required=True)
    ap.add_argument("--output-csv",type=Path,required=True)
    ap.add_argument("--output-json",type=Path,required=True)
    args=ap.parse_args()

    contract=load_contract(args.contract)
    if sha256_path(args.survivor_geometry_csv)!=EXPECTED_SURVIVOR_GEOMETRY_SHA256:
        raise RuntimeError("Phase-2 survivor geometry CSV SHA256 drift")
    if sha256_path(args.authorization_json)!=EXPECTED_AUTHORIZATION_SHA256:
        raise RuntimeError("Phase-4 authorization SHA256 drift")

    source=json.loads(args.realm_source_receipt.read_text())
    if source.get("schema")!="ttf_palearctic_realm_source_receipt_v0.1":
        raise RuntimeError("unexpected realm-source receipt")
    if source.get("status")!="FROZEN_BEFORE_PANEL_MEMBERSHIP":
        raise RuntimeError("realm source was not frozen before panel membership")
    subset=json.loads(args.realm_subset_receipt.read_text())
    if subset.get("schema")!="ttf_palearctic_realm_subset_receipt_v0.1":
        raise RuntimeError("unexpected Palearctic realm subset receipt")
    if subset.get("status")!="DERIVED_BEFORE_PANEL_MEMBERSHIP":
        raise RuntimeError("Palearctic subset was not derived before panel membership")
    if subset.get("source_receipt_sha256")!=sha256_path(args.realm_source_receipt):
        raise RuntimeError("Palearctic subset/source receipt linkage drift")
    if sha256_path(args.ecoregions)!=subset["geojson_sha256"]:
        raise RuntimeError("Palearctic GeoJSON SHA256 differs from frozen subset receipt")

    auth=json.loads(args.authorization_json.read_text())
    survivors=set(auth["species"]["train_species"]+auth["species"]["eval_species"])
    if len(survivors)!=211:
        raise RuntimeError("authorization survivor count drift")

    with args.survivor_geometry_csv.open(newline="",encoding="utf-8") as h:
        rows=list(csv.DictReader(h))
    required={"species","class","order","family","latitude","longitude"}
    if not rows or not required.issubset(rows[0]):
        raise RuntimeError("Phase-2 survivor geometry CSV schema drift")

    by={}
    for row in rows:
        sp=str(row["species"])
        if sp not in survivors:
            continue
        by.setdefault(sp,[]).append(row)
    if set(by)!=survivors:
        missing=sorted(survivors-set(by))
        raise RuntimeError(f"Phase-2 survivor geometry is missing authorized survivors: {missing[:10]}")

    try:
        import geopandas as gpd
        from shapely.geometry import Point
        from shapely.ops import unary_union
    except ImportError as exc:
        raise RuntimeError("geopandas and shapely are required for the realm screen") from exc

    eco=gpd.read_file(args.ecoregions)
    realm_field=contract["panel_rule"]["realm_source"]["realm_field"]
    realm_value=contract["panel_rule"]["realm_source"]["realm_value"]
    if realm_field not in eco.columns:
        raise RuntimeError(f"ecoregions layer lacks required realm field {realm_field}")
    if eco.crs is None:
        raise RuntimeError("ecoregions layer lacks CRS")
    eco=eco.to_crs("EPSG:4326")
    pa=eco.loc[eco[realm_field].astype(str)==realm_value]
    if pa.empty:
        raise RuntimeError("Palearctic realm selection returned zero polygons")
    pa_union=unary_union(list(pa.geometry))

    allowed=set(contract["panel_rule"]["taxonomy"]["allowed_orders"])
    threshold=float(
        contract["panel_rule"]["palearctic_membership"]["minimum_fraction_inside_palearctic"]
    )

    output=[]
    for sp in sorted(survivors):
        sp_rows=by[sp]
        cls=str(sp_rows[0]["class"])
        order=str(sp_rows[0]["order"])
        family=str(sp_rows[0]["family"])
        coords=[]
        for row in sp_rows:
            try:
                lat=float(row["latitude"]); lon=float(row["longitude"])
            except ValueError:
                coords.append((None,None,False))
                continue
            inside=bool(pa_union.covers(Point(lon,lat)))
            coords.append((lat,lon,inside))
        n=len(coords)
        inside=sum(1 for _,_,yes in coords if yes)
        fraction=inside/n if n else 0.0
        taxon_ok=cls=="Insecta" and order in allowed
        eligible=taxon_ok and fraction>=threshold
        output.append({
            "species":sp,
            "class":cls,
            "order":order,
            "family":family,
            "phase1_localities":n,
            "palearctic_localities":inside,
            "palearctic_fraction":fraction,
            "taxonomy_eligible":taxon_ok,
            "eligible":eligible,
        })

    eligible=sorted(row["species"] for row in output if row["eligible"])
    minimum=int(contract["panel_rule"]["minimum_species_to_continue"])
    source_species,target_species=deterministic_split(
        eligible,
        "palearctic-holometabola-lgm-v0.1",
    )
    status=(
        "PASS_PANEL_SIZE_TO_LGM_SDM"
        if len(eligible)>=minimum
        else "STOP_NOT_EVALUABLE_FOR_LGM_PAIRWISE_PROGRAM"
    )

    args.output_csv.parent.mkdir(parents=True,exist_ok=True)
    with args.output_csv.open("w",newline="",encoding="utf-8") as h:
        writer=csv.DictWriter(h,fieldnames=list(output[0]))
        writer.writeheader(); writer.writerows(output)

    payload={
        "schema":"ttf_genetic_palearctic_holometabola_panel_v0.1",
        "status":status,
        "contract_sha256":sha256_path(args.contract),
        "survivor_geometry_csv_sha256":sha256_path(args.survivor_geometry_csv),
        "authorization_sha256":sha256_path(args.authorization_json),
        "realm_source_receipt_sha256":sha256_path(args.realm_source_receipt),
        "realm_subset_receipt_sha256":sha256_path(args.realm_subset_receipt),
        "palearctic_geojson_sha256":sha256_path(args.ecoregions),
        "survivor_species":len(survivors),
        "eligible_species":len(eligible),
        "minimum_species_to_continue":minimum,
        "source_species":source_species,
        "target_species":target_species,
        "source_count":len(source_species),
        "target_count":len(target_species),
        "primary_dyads":len(source_species)*len(target_species),
        "genetic_pair_response_opened":False,
        "species_rows":output,
    }
    args.output_json.parent.mkdir(parents=True,exist_ok=True)
    args.output_json.write_text(json.dumps(payload,indent=2,sort_keys=True)+"\n")
    print(json.dumps({
        "status":status,
        "eligible_species":len(eligible),
        "source_count":len(source_species),
        "target_count":len(target_species),
        "primary_dyads":len(source_species)*len(target_species),
        "genetic_pair_response_opened":False,
    },sort_keys=True))
    return 0


if __name__=="__main__":
    raise SystemExit(main())
