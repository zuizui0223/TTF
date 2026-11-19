"""Cross-paper co-location QA: Raabova 2007 population sites via Sudova 2010 Table 1.

Source-only. The published DMS site coordinates are precommitted. Hargreaves
candidate source positions are original same-study/species GARDEN site
coordinates where exactly matching source and garden population IDs exist.
Nearest site distance validates spatial concordance, not primary experiment ID.
No biological fitness or historical climate outcome opened or exported.
"""
from __future__ import annotations
import argparse,hashlib,json,re
from pathlib import Path
from bontrager_hargreaves_finegrain_overlap_v01 import original_bytes, haversine
from hargreaves_source_as_garden_georecovery_v01 import source_metadata, one_point

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
CONTRACT=HERE/"raabova_sudova_crosspaper_sites_contract_v01.json"
BASE_CONTRACT=HERE/"hargreaves_source_as_garden_georecovery_contract_v01.json"

def dms(s:str, axis:str)->float:
    m=re.fullmatch(r"(\d{1,3}):(\d{1,2}):(\d{1,2}(?:\.\d+)?)([NSEW])",s)
    if not m:
        raise ValueError("Not an original verified DMS coordinate")
    deg,minutes,seconds,compass=m.groups()
    deg=int(deg);minutes=int(minutes);seconds=float(seconds)
    if minutes>=60 or seconds>=60:
        raise ValueError("Invalid DMS minute or second")
    if axis=="lat":
        if compass not in ("N","S") or deg>90:
            raise ValueError("Latitude hemisphere/range drift")
    elif axis=="lon":
        if compass not in ("E","W") or deg>180:
            raise ValueError("Longitude hemisphere/range drift")
    else:
        raise ValueError("Unknown axis")
    val=deg+minutes/60.+seconds/3600.
    return -val if compass in ("S","W") else val

def primary_sites(contract:dict):
    places=contract["original_external_localities"]
    if len(places)!=4 or [x["source_site_locality_number_in_Raabova_2007"] for x in places]!=[3,1,5,6]:
        raise ValueError("Frozen external primary source localities/order changed")
    return [(
        v["source_site_locality_number_in_Raabova_2007"],
        (dms(v["longitude_dms"],"lon"),dms(v["latitude_dms"],"lat"))
    ) for v in places]

def colocate(primary:list, candidates:list,contract:dict)->dict:
    rule=contract["prespecified_distance_rule"]
    limit=rule["maximum_distance_for_cross_paper_corroboration_km"]
    margin=rule["minimum_nearest_second_best_separation_km"]
    if (limit!=1.0 or margin!=1.0
        or not rule["each_hargreaves_site_can_match_at_most_one_external_site"]):
        raise ValueError("No post-source position tolerance tuning")
    if len(set(candidates))!=len(candidates):
        raise ValueError("Nonunique original source-as-garden coordinates")
    results=[];used=set()
    for siteid,original_xy in primary:
        ds=sorted((haversine(original_xy,xy),i) for i,xy in enumerate(candidates))
        if not ds:
            raise ValueError("No same-study source-as-garden sites")
        d1,ix=ds[0]
        d2=ds[1][0] if len(ds)>1 else None
        unique=(d2 is None or d2-d1>=margin)
        qualifying=d1<=limit and unique and ix not in used
        if qualifying:used.add(ix)
        bucket=("<0.1km" if d1<0.1 else "0.1-1km" if d1<=1 else ">1km")
        results.append({
            "original_primary_locality_number":siteid,
            "source_only_nearest_distance_bucket":bucket,
            "unique_nearest_population_1km_corroboration":qualifying,
            "nearest_distance_km":round(d1,5),
            "nearest_second_choice_ambiguous":not unique,
            "actual_Hargreaves_source_ID_or_coordinates_exported":False
        })
    return {
        "original_external_locality_count":len(primary),
        "Hargreaves_original_same_study_source_as_garden_population_positions":len(candidates),
        "published_sites_spatially_corroborated_uniquely_within_1km":sum(x["unique_nearest_population_1km_corroboration"] for x in results),
        "published_sites_ambiguous_or_unmatched":sum(not x["unique_nearest_population_1km_corroboration"] for x in results),
        "primary_paper_exact_population_ID_to_Hargreaves_source_ID_crosswalk_certified":False,
        "original_field_GPS_independent_acquisition_certified":False,
        "external_original_site_match_details":results,
    }

def run():
    p=argparse.ArgumentParser()
    p.add_argument("--output",type=Path,required=True)
    args=p.parse_args()
    contract=json.loads(CONTRACT.read_text(encoding="utf-8"))
    if contract["status"]!="FROZEN_BEFORE_HARGREAVES_SOURCE_AS_GARDEN_CROSSPAPER_SITE_CHECK":
        raise ValueError("No frozen external original literature site-authority")
    prior=ROOT/contract["source_as_garden_original"]["prior_source_as_garden_receipt"]
    rawprior=prior.read_bytes()
    if hashlib.sha256(rawprior).hexdigest()!=contract["source_as_garden_original"]["prior_receipt_SHA256"]:
        raise ValueError("Original source-as-garden result receipt differs")
    priorinfo=json.loads(rawprior)
    if (not priorinfo["no_lgm_climate_or_genetic_fitness_effect_opened"]
        or not priorinfo["Bontrager_frozen_44_of_50_gate_unchanged"]):
        raise ValueError("Original source-only project scope drift")
    source_rule=json.loads(BASE_CONTRACT.read_text())
    main=json.loads((ROOT/source_rule["original_source_contract"]).read_text())
    original=original_bytes(main["original"]["hargreaves"])
    if hashlib.sha256(original).hexdigest()!=contract["source_as_garden_original"]["exact_source_SHA256"]:
        raise ValueError("Original Hargreaves source hash differs")
    obs=source_metadata(original,"hargreaves",source_rule)
    study=contract["target_study_label"];taxon=contract["taxon"]
    positions={one_point(obs["gardens"].get(k,set())) for k in obs["origin_keys"]
               if k[0]==study and k[1]==taxon and k in obs["gardens"]}
    positions.discard(None)
    if len(positions)!=contract["source_as_garden_original"]["expected_study_recovered_unique_source_sites"]:
        raise ValueError("Study species source-as-garden input geometry drift")
    result=colocate(primary_sites(contract),sorted(positions),contract)
    out={
        "schema":"ttf_Raabova_Sudova_pub_primary_sites_crosspaper_source_only_geocoherence_v01",
        "status":"LITERATURE_SAME_SITE_SPATIAL_COHERENCE_NOT_IDENTICAL_SOURCE_ID_CERTIFICATION",
        "original_study_label":study,
        "original_primary_study_doi":contract["original_transplant_study_doi"],
        "cross_paper_site_source_doi":contract["external_primary_source_doi"],
        "source_as_garden_receipt_SHA256":contract["source_as_garden_original"]["prior_receipt_SHA256"],
        "source_only_primary_published_site_correspondence":result,
        "original_population_lon_lat_values_exported":False,
        "raw_original_source_pop_ID_exported":False,
        "fitness_or_adaptation_effect_used":False,
        "historic_climate_measurement_fitted":False,
        "Bontrager_frozen_44_of_50_gate_unchanged":True
    }
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(out,indent=2,ensure_ascii=False,sort_keys=True)+"\n")
    print("RAABOVA_PRIMARY_CROSSPAPER_GEO_VALIDATION=",json.dumps({
        "status":out["status"],"source_study":study,
        "summary":{k:v for k,v in result.items() if k!="external_original_site_match_details"}
    },sort_keys=True))

if __name__=="__main__":
    run()
