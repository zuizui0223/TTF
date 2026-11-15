"""Population-as-garden coordinate recovery: strictly source-metadata only.

Discover whether a reciprocal-transplant origin ID also occurs as an EXACT
same-study/same-species garden ID with one unique valid destination coordinate.
Benchmark those candidates on previously known Bontrager original source
coordinates; NEVER use response values or fit historical climate.
"""
from __future__ import annotations
import argparse,csv,hashlib,io,json,statistics
from collections import Counter,defaultdict
from pathlib import Path
from bontrager_hargreaves_finegrain_overlap_v01 import (
    original_bytes, read_bontrager_csv, norm_taxon, valid_coord, haversine
)

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
CONTRACT=HERE/"hargreaves_source_as_garden_georecovery_contract_v01.json"

def source_metadata(raw:bytes,cohort:str, rule:dict) -> dict:
    fields=rule["source_only_allowlist_"+cohort]
    required_b=["study","functional_group","taxon","source","site","lat_source","long_source","lat_site","long_site"]
    required_h=["study","functional_group","taxon","source","site","experiment","fitness_type","lat_site","long_site"]
    if fields!=(required_b if cohort=="bontrager" else required_h):
        raise ValueError("Source-only original metadata field allowlist drifted")
    if cohort=="hargreaves":
        decoded=raw.decode("cp1252",errors="strict")
    elif cohort=="bontrager":
        decoded=raw.decode("utf-8-sig",errors="strict")
    else:
        raise ValueError("Unknown original source")
    reader=csv.reader(io.StringIO(decoded,newline=""), strict=True)
    header=next(reader,[])
    if len(header)!=47 or len(set(header))!=47 or not set(fields).issubset(header):
        raise ValueError("Original CSV 47-column metadata source changed")
    at={key:header.index(key) for key in fields}
    plant_groups=set(rule["plant_functional_groups"])
    if plant_groups!={"annual","herbaceous_perennial","woody_perennial"}:
        raise ValueError("Original plant functional groups drifted")
    origin_keys=set()
    gardens=defaultdict(set)
    original_sources=defaultdict(set)
    cells=defaultdict(set)
    nrows=0;plant_rows=0
    for row in reader:
        if not row:continue
        nrows+=1
        if nrows>200000 or len(row)!=47:
            raise ValueError("Original row count/width changed")
        # Only fields declared above. NO fitness response access.
        m={f:row[pos].strip() for f,pos in at.items()}
        if m["functional_group"] not in plant_groups:continue
        plant_rows+=1
        study=m["study"]
        taxon=norm_taxon(m["taxon"])
        if not study or not taxon:continue
        base=(study,taxon)
        s=m["source"];g=m["site"]
        if s:origin_keys.add(base+(s,))
        if g:
            gardens[base+(g,)].add(valid_coord(m["long_site"],m["lat_site"]))
        if cohort=="bontrager" and s:
            original_sources[base+(s,)].add(valid_coord(m["long_source"],m["lat_source"]))
        if cohort=="hargreaves" and s and g and m["fitness_type"]:
            exp=m["experiment"] or "<UNKNOWN_EXPERIMENT>"
            cells[base+(exp,g,m["fitness_type"])].add(s)
    return {"cohort":cohort,"original_rows":nrows,"plant_rows":plant_rows,
            "origin_keys":origin_keys,"gardens":gardens,
            "source_truth":original_sources,"source_groups":cells}

def one_point(values:set)->tuple|None:
    # Treat either coordinate inconsistency OR ANY missing marker as ambiguous.
    if len(values)==1 and None not in values:
        return next(iter(values))
    return None

def percent_at(data, q):
    if not data:return None
    a=sorted(data)
    i=(len(a)-1)*q
    lo=int(i);hi=min(lo+1,len(a)-1)
    return a[lo]+(a[hi]-a[lo])*(i-lo)

def evaluation(h:dict,b:dict,rule:dict,receipt:dict) -> dict:
    frozen=rule["study_label_partition"]
    targets=receipt["Hargreaves_unlinked_full_study_labels_for_followup"]
    if sorted(targets)!=sorted(frozen.get("exact_unlinked_labels",targets)):
        raise ValueError("Original unlinked source study labels changed")
    if len(targets)!=frozen["study_count"] or len(set(targets))!=35:
        raise ValueError("Prior frozen Hargreaves-only 35 study labels not reproduced")
    known_source_positions={k:one_point(v) for k,v in b["source_truth"].items()}
    recovered={k:one_point(h["gardens"].get(k,set())) for k in h["origin_keys"]}
    shared=set(h["origin_keys"])&set(known_source_positions)
    source_as_garden_keys={k for k in h["origin_keys"] if k in h["gardens"]}
    valid_keys={k for k in source_as_garden_keys if recovered[k] is not None}
    distances=[]
    calibration_totals=Counter()
    for k in sorted(shared):
        truth=known_source_positions[k]
        predicted=recovered.get(k)
        if predicted is None:
            calibration_totals["no_unique_eligible_garden_coordinate"]+=1
            continue
        if truth is None:
            calibration_totals["bontrager_original_truth_missing_or_ambiguous"]+=1
            continue
        calibration_totals["comparable_source_origin_keys"]+=1
        distance=haversine(predicted,truth)
        distances.append(distance)
        if distance<=1:calibration_totals["km_le_1"]+=1
        elif distance<=10:calibration_totals["km_gt_1_le_10"]+=1
        else:calibration_totals["km_gt_10"]+=1
    per_study=[]
    for study in sorted(targets):
        keys={k for k in h["origin_keys"] if k[0]==study}
        direct={k for k in keys if k in source_as_garden_keys}
        good={k for k in direct if k in valid_keys}
        geometry=[(k,src) for k,src in h["source_groups"].items() if k[0]==study and len(src)>=2]
        recoverable=[]
        for cell, sources in geometry:
            key_prefix=cell[:2]
            geo={recovered[key_prefix+(s,)] for s in sources
                 if recovered.get(key_prefix+(s,)) is not None}
            if len(geo)>=2:
                recoverable.append(cell)
        per_study.append({
            "original_study_label":study,
            "distinct_source_population_ID_keys":len(keys),
            "source_ID_exactly_present_as_same_study_species_garden_ID":len(direct),
            "with_unique_valid_original_destination_coordinate":len(good),
            "unique_georecovered_source_as_garden_coordinates":len({recovered[k] for k in good}),
            "source_only_original_multi_source_design_cells":len(geometry),
            "multi_source_cells_with_two_distinct_georecovered_coordinates":len(recoverable),
            "original_source_geocoordinates_independently_verified":False,
            "primary_original_experiment_independence_certified":False
        })
    summary={
        "original_Hargreaves_source_population_ID_keys_all_plant":len(h["origin_keys"]),
        "original_Hargreaves_source_ID_keys_existing_as_same_species_gardens":len(source_as_garden_keys),
        "original_Hargreaves_source_ID_keys_with_unique_garden_coordinate":len(valid_keys),
        "shared_Hargreaves_Bontrager_exact_source_population_keys":len(shared),
        "calibration_true_source_comparisons":len(distances),
        "calibration_distance_median_km":percent_at(distances,0.5),
        "calibration_distance_p95_km":percent_at(distances,0.95),
        "calibration_distance_max_km":max(distances) if distances else None,
        "calibration_distance_bins":dict(sorted(calibration_totals.items())),
        "calibration_minimum_key_count_for_reporting_met":
            len(distances)>=rule["strict_source_as_garden_rule"]["minimum_calibration_population_keys_to_make_any_reliability_statement"],
        "unlinked_Hargreaves_original_study_labels":len(targets),
        "unlinked_study_labels_with_any_pop_as_garden_recovered":sum(x["with_unique_valid_original_destination_coordinate"]>0 for x in per_study),
        "unlinked_study_labels_with_two_or_more_pop_as_garden_recovered":sum(x["unique_georecovered_source_as_garden_coordinates"]>=2 for x in per_study),
        "unlinked_study_labels_with_any_same_garden_two_georecovered_origins":
            sum(x["multi_source_cells_with_two_distinct_georecovered_coordinates"]>0 for x in per_study),
        "unlinked_design_cells_with_two_georecovered_origins":
            sum(x["multi_source_cells_with_two_distinct_georecovered_coordinates"] for x in per_study),
        "unlinked_distinct_source_IDs_eligible_for_garden_derived_proxy":
            sum(x["with_unique_valid_original_destination_coordinate"] for x in per_study),
        "independently_verified_original_source_latlon_from_original_Hargreaves_columns":0,
        "independently_certified_original_source_experiments":0
    }
    return {
        "schema":"ttf_Hargreaves_reciprocal_exact_population_as_garden_georecovery_source_only_v01",
        "status":"SOURCE_AS_GARDEN_ONLY_DESIGN_CANDIDATES_NOT_ORIGINAL_GEO_OR_ADAPTATION",
        "original_source_row_counts":{"hargreaves":h["original_rows"],"bontrager":b["original_rows"]},
        "summary":summary,"per_unlinked_study":per_study,
        "output_contains_raw_coordinates_or_source_IDs":False,
        "original_garden_coordinate_proxy_is_NOT_original_source_coordinate_certification":True,
        "no_lgm_climate_or_genetic_fitness_effect_opened":True,
        "Bontrager_frozen_44_of_50_gate_unchanged":True
    }

def run():
    p=argparse.ArgumentParser()
    p.add_argument("--output",required=True,type=Path)
    args=p.parse_args()
    rule=json.loads(CONTRACT.read_text())
    if rule["status"]!="FROZEN_BEFORE_NEW_SOURCE_AS_GARDEN_LABEL_JOIN":
        raise ValueError("Source-only reciprocal-georecovery rule not frozen")
    original_contract=json.loads((ROOT/rule["original_source_contract"]).read_text())
    rpath=ROOT/rule["study_label_partition"]["previously_unlinked_35_receipt"]
    rbytes=rpath.read_bytes()
    if hashlib.sha256(rbytes).hexdigest()!=rule["study_label_partition"]["previously_unlinked_35_receipt_sha256"]:
        raise ValueError("Original unlinked 35 receipt SHA differs")
    receipt=json.loads(rbytes)
    x=original_contract["original"]
    bontrager=original_bytes(x["bontrager"])
    hargreaves=original_bytes(x["hargreaves"])
    b=source_metadata(read_bontrager_csv(bontrager,x["bontrager"]["member"]),"bontrager",rule)
    h=source_metadata(hargreaves,"hargreaves",rule)
    if (b["original_rows"],h["original_rows"])!=(7544,11088):
        raise ValueError("Previously verified original source physical row counts differ")
    result=evaluation(h,b,rule,receipt)
    result["authentic_original_source_sha256"]={
        "bontrager":x["bontrager"]["sha256"],"hargreaves":x["hargreaves"]["sha256"]}
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(result,indent=2,ensure_ascii=False,sort_keys=True)+"\n",encoding="utf-8")
    print("SOURCE_AS_GARDEN_RECIPROCAL_GEO=",json.dumps(result["summary"],sort_keys=True))

if __name__=="__main__":
    run()
