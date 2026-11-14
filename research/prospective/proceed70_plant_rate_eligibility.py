"""Source-only trait-family/time-denominator inventory for PROCEED70 plants.

All biological-rate/mean/variance fields remain unread. Counts use the same
fixed 70 georeferenced species and disjoint 35/35 study-species components as
the original 253-site CHELSA8 source. Does not choose a phenotype outcome.
"""
from __future__ import annotations

import argparse,csv,hashlib,io,json
from collections import Counter,defaultdict
from pathlib import Path

from proceed_v62_metadata_row_census import (
    fetch_original_datafile_bytes,valid_binomial,valid_coordinates)
from proceed_plant_chelsa_source_support import source_geometry,component_assignment
from proceed_v62_genetic_geometry_intersection import ROW_CONTRACT,SOURCE_CONTRACT,NATIVE_LOCK
from proceed_v62_borealis_header_only import validated_get_catalog,verify_frozen_catalog

HERE=Path(__file__).resolve().parent
CONTRACT=HERE/"proceed70_plant_rate_eligibility_contract_v01.json"


def safe_count_labels(z,max_codes=30):
    if len(z)>max_codes:
        return {"distinct_source_codes":len(z),"too_many_codes_to_report":True}
    return dict(sorted((str(k or "<missing>"),int(v)) for k,v in z.items()))


def response_blind_rate_metadata(raw,contract,check_authoritative=True):
    allowed=contract["metadata_fields_allowed"]
    prohibited=set(contract["explicitly_forbidden_response_fields"])
    if set(allowed)&prohibited:
        raise ValueError("Forbidden rate data in source metadata allowlist")
    r=csv.reader(io.TextIOWrapper(io.BytesIO(raw),encoding="utf-8-sig",newline=""),
                 delimiter=",",strict=True)
    header=next(r)
    if len(header)!=69 or len(set(header))!=69 or not set(allowed).issubset(header):
        raise ValueError("Original v6.2 69-field metadata source schema changed")
    ix={k:header.index(k) for k in allowed}
    original_rows=0;metarows=0;study_species=defaultdict(set)
    families=defaultdict(lambda: {
        "source_rows":0,"species":set(),"studies":set(),"genera":set(),
        "geography":set(),"design":Counter(),"datatype":Counter(),
        "experimental":Counter(),"exp_type":Counter(),
        "years_present":0,"generations_present":0,"gtime_present":0,
    })
    overall_design=Counter();overall_datatype=Counter();overall_experimental=Counter()
    original_genetics=Counter();all_trait_types=Counter()
    for values in r:
        if not values:continue
        original_rows+=1
        if len(values)!=69:
            raise ValueError("Source original native CSV record width changed")
        m={k:values[i].strip() for k,i in ix.items()}
        if m["phylum"]!="Streptophyta" or m["taxa"]!="Plant" or m["genphen"]!="Genetic":
            continue
        c1=valid_coordinates(m["sample1_latitude"],m["sample1_longitude"])
        c2=valid_coordinates(m["sample2_latitude"],m["sample2_longitude"])
        sp=valid_binomial(m["sp_pub"]) or valid_binomial(m["sp_ncbi"])
        study=m["ref_ID"]
        if not(c1 and c2 and sp and study):
            continue
        metarows+=1
        study_species[study].add(sp)
        code=m["trait_type"] or "<missing>"
        a=families[code]
        a["source_rows"]+=1
        a["species"].add(sp);a["studies"].add(study);a["genera"].add(sp.split()[0])
        # Original study/species/ordered locality pair, independent of
        # how many phenotypic traits/effect rows used that geography.
        a["geography"].add((study,sp,c1,c2))
        a["design"][m["design"]]+=1
        a["datatype"][m["data_type"]]+=1
        a["experimental"][m["experimental"]]+=1
        a["exp_type"][m["exp_type"]]+=1
        a["years_present"]+=bool(m["years"])
        a["generations_present"]+=bool(m["generations"])
        a["gtime_present"]+=bool(m["gtime"])
        all_trait_types[code]+=1
        overall_design[m["design"]]+=1
        overall_datatype[m["data_type"]]+=1
        overall_experimental[m["experimental"]]+=1

    roles,split,digest=component_assignment(study_species)
    expected=contract["authoritative_inputs"]
    if check_authoritative and (original_rows!=9263 or metarows!=expected["georeferenced_eligible_source_rows"]
        or len(roles)!=70 or len(study_species)!=104
        or digest!=expected["study_species_role_sha256"]):
        raise ValueError("Frozen original plant source population or role SHA mismatch")
    result={}
    for code,v in sorted(families.items()):
        groups={k:len({s for s in v["species"] if roles[s]==k})
                for k in ("train","heldout")}
        result[code]={
            "original_rate_rows_counted_as_metadata_only":v["source_rows"],
            "distinct_original_species":len(v["species"]),
            "distinct_original_studies":len(v["studies"]),
            "distinct_original_genera":len(v["genera"]),
            "distinct_study_species_ordered_geo_pairs":len(v["geography"]),
            "species_by_original_study_species_disjoint_role":groups,
            "source_design_codes":safe_count_labels(v["design"]),
            "data_type_codes":safe_count_labels(v["datatype"]),
            "source_experimental_codes":safe_count_labels(v["experimental"]),
            "experimental_type_codes":safe_count_labels(v["exp_type"]),
            "years_metadata_present":v["years_present"],
            "generations_metadata_present":v["generations_present"],
            "generation_time_metadata_present":v["gtime_present"],
        }
    return {
        "schema":"ttf_proceed70_response_blind_trait_family_and_rate_design_source_census_v01",
        "status":"TRAIT_AND_TIME_SOURCE_CATEGORIES_CENSUSED_NO_RATE_OUTCOME_OPENED",
        "source_total_original_records":original_rows,
        "original_eligible_plant_effect_source_records":metarows,
        "original_study_species_role_sha256":digest,
        "original_train_species":split["species_per_role"]["train"],
        "original_heldout_species":split["species_per_role"]["heldout"],
        "trait_family_metadata":result,
        "source_data_type_global_counts":safe_count_labels(overall_datatype),
        "original_study_design_global_counts":safe_count_labels(overall_design),
        "original_experimental_global_counts":safe_count_labels(overall_experimental),
        "phenotypic_change_or_genetic_effect_rate_values_read":False,
        "model_fit_or_history_effect_tested":False,
        "original_TTF_outcomes_reopened":False,
        "original_source_coordinate_or_trait_rows_exported":False,
        "source_genetic_label_does_not_establish_adaptive_fitness":True,
    }


def run():
    parser=argparse.ArgumentParser()
    parser.add_argument("--output",type=Path,required=True)
    args=parser.parse_args()
    c=json.loads(CONTRACT.read_text())
    if c["status"]!="FROZEN_SOURCE_ONLY_TRAIT_AND_TEMPORAL_DESIGN_ELIGIBILITY_BEFORE_ANY_CHANGE_RATE":
        raise ValueError("Source-only trait rate qualification plan not frozen")
    src=json.loads(SOURCE_CONTRACT.read_text())["current_source"]
    native=json.loads(NATIVE_LOCK.read_text())["native_original_transport"]
    out={"schema":"ttf_proceed70_response_blind_trait_family_and_rate_design_source_census_v01",
         "phenotypic_change_or_genetic_effect_rate_values_read":False,
         "model_fit_or_history_effect_tested":False,
         "original_TTF_outcomes_reopened":False,
         "original_source_coordinate_or_trait_rows_exported":False}
    try:
        cat=validated_get_catalog("https://borealisdata.ca/api/datasets/:persistentId/?persistentId=doi%3A10.5683%2FSP3%2FNXSL3Q")
        verify_frozen_catalog(cat,src)
        raw,md5,sha=fetch_original_datafile_bytes(1024206)
        if (md5!=native["observed_md5"] or sha!=native["observed_sha256"]
            or len(raw)!=native["observed_size_bytes"]):
            raise ValueError("Original source rate-row bytes changed")
        out.update(response_blind_rate_metadata(raw,c))
        out["source_sha256"]=sha
        out["source_md5"]=md5
    except Exception as exc:
        out["status"]="NOT_EVALUABLE_TRAIT_FAMILY_METADATA_SOURCE_OR_REPLAY"
        out["error_type"]=type(exc).__name__
        out["error_message"]=str(exc)[:400]
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(out,indent=2,sort_keys=True)+"\n")
    print(json.dumps({k:out.get(k) for k in (
        "status","original_eligible_plant_effect_source_records",
        "source_data_type_global_counts","original_study_design_global_counts")},
        sort_keys=True))


if __name__=="__main__":
    run()
