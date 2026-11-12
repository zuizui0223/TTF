"""Source-only plant genetic-category study/species/coordinate macro qualification.

Reads original SHA256/MD5-bound PROCEED v6.2 native CSV, but accesses ONLY
predeclared metadata columns. Genetic trait rates, effect-size values, means,
standard deviations and actual genotypes are never interpreted or fitted.
"""
from __future__ import annotations
import argparse,csv,hashlib,io,json
from collections import Counter,defaultdict
from pathlib import Path

from proceed_v62_borealis_header_only import validated_get_catalog,verify_frozen_catalog
from proceed_v62_metadata_row_census import (
    fetch_original_datafile_bytes,valid_binomial,valid_coordinates,safe_counts)
from proceed_v62_genetic_geometry_intersection import SOURCE_CONTRACT,NATIVE_LOCK,ROW_CONTRACT

HERE=Path(__file__).resolve().parent
CONTRACT=HERE/"proceed_v62_plant_genetic_geo_contract_v01.json"


def component_member_sets(study_species):
    species_study=defaultdict(set)
    for study,sps in study_species.items():
        for s in sps:species_study[s].add(study)
    pending=set(species_study)
    results=[]
    while pending:
        todo=[min(pending)]; ss=set();rr=set()
        while todo:
            s=todo.pop()
            if s in ss:continue
            ss.add(s);pending.discard(s)
            for study in sorted(species_study[s]):
                if study in rr:continue
                rr.add(study)
                todo.extend(study_species[study]-ss)
        results.append((ss,rr))
    return results


def split_connected_components(study_species):
    """Response-blind deterministic component-blocked split; studies and species disjoint."""
    comps=component_member_sets(study_species)
    def key(comp):
        name="|".join(sorted(comp[0]))
        h=hashlib.sha256(("ttf-proceed-v62-plant-component-v01|"+name).encode()).hexdigest()
        return (-len(comp[0]),h)
    comps.sort(key=key)
    role={"train":set(),"heldout":set()}
    role_study={"train":set(),"heldout":set()}
    for taxa,studies in comps:
        role_key="train" if len(role["train"])<=len(role["heldout"]) else "heldout"
        role[role_key].update(taxa)
        role_study[role_key].update(studies)
    if role["train"]&role["heldout"] or role_study["train"]&role_study["heldout"]:
        raise AssertionError("Study/species graph leakage")
    component_sizes=[{"species":len(s),"studies":len(t)} for s,t in comps]
    return {
        "components":len(comps),
        "component_sizes_first_12":component_sizes[:12],
        "species_per_role":{k:len(v) for k,v in role.items()},
        "studies_per_role":{k:len(v) for k,v in role_study.items()},
        "roles_study_species_disjoint":True,
        "source_only_component_membership_sha256":hashlib.sha256(
            json.dumps({k:sorted(v) for k,v in role.items()},sort_keys=True,separators=(",",":")).encode()
        ).hexdigest(),
        "exposes_species_list":False
    }


def source_only_plant_counts(raw,contract):
    permitted=json.loads(ROW_CONTRACT.read_text())["source_metadata_columns_authorized"]
    forbidden=set(json.loads(ROW_CONTRACT.read_text())["outcome_fields_must_not_be_decoded"])
    if set(permitted)&forbidden:raise ValueError("Data response column in metadata allowlist")
    reader=csv.reader(io.TextIOWrapper(io.BytesIO(raw),encoding="utf-8-sig",newline=""),delimiter=",",strict=True)
    header=next(reader)
    if len(header)!=69 or len(set(header))!=69 or not set(permitted).issubset(header):
        raise ValueError("Original native CSV header mismatch")
    ix={k:header.index(k) for k in permitted}
    total=0;plantgen=0;geocoded=0;missingids=0;invalidspecies=0
    original_studies=set();plants=set();genera=set();systems=set()
    study_species=defaultdict(set)
    design=Counter();experimental=Counter();exp_type=Counter();trait=Counter();time_label=Counter()
    coordinate_type=Counter();source_distinct_coord=set()
    group_duplicate=Counter()
    for fields in reader:
        if not fields:continue
        total+=1
        if len(fields)!=69:raise ValueError("Original source row width not reproducible")
        # Do not access a single forbidden biological-response field.
        m={name:fields[i].strip() for name,i in ix.items()}
        if m["phylum"]!=contract["source_eligibility"]["phylum"] or m["genphen"]!=contract["source_eligibility"]["genphen"]:
            continue
        plantgen+=1
        a=valid_coordinates(m["sample1_latitude"],m["sample1_longitude"])
        b=valid_coordinates(m["sample2_latitude"],m["sample2_longitude"])
        coordinate_type["both" if a and b else "first_only" if a else "second_only" if b else "neither"]+=1
        if not (a and b):continue
        geocoded+=1
        name=valid_binomial(m["sp_pub"]) or valid_binomial(m["sp_ncbi"])
        if not name:
            invalidspecies+=1;continue
        study=m["ref_ID"]
        if not study:
            missingids+=1;continue
        if m["taxa"]!=contract["source_eligibility"]["taxa_group_frozen"]:
            raise ValueError("Plant phylum/genetic source has unexpected taxonomic group label")
        plants.add(name);genera.add(name.split()[0]);original_studies.add(study)
        study_species[study].add(name)
        systems.add(m["sys_ID"])
        group_duplicate[(study,name)]+=1
        design[m["design"]]+=1
        experimental[m["experimental"]]+=1
        exp_type[m["exp_type"]]+=1
        trait[m["trait_type"]]+=1
        time_label["years_present" if m["years"] else "years_missing"]+=1
        time_label["generations_present" if m["generations"] else "generations_missing"]+=1
        source_distinct_coord.add((a,b))
    split=split_connected_components(study_species)
    constraints=contract["source_eligibility"]
    macro=(len(plants)>=constraints["minimum_independent_taxa"]
           and len(original_studies)>=constraints["minimum_source_studies"]
           and len(genera)>=constraints["minimum_genera"]
           and min(split["species_per_role"].values())>=constraints["minimum_species_in_each_study_species_disjoint_role"])
    return {
      "schema":"ttf_proceed_v62_original_plant_genetic_change_source_eligibility_v01",
      "status":"SOURCE_ONLY_MACRO_GEOMETRY_POSSIBLE_NOT_ADAPTATION" if macro else "PLANT_GENETIC_GEO_SOURCE_NOT_MACRO_EVALUABLE",
      "all_original_source_rows":total,
      "plant_genetic_category_rows_all_geographic_status":plantgen,
      "plant_genetic_geolocation_status":dict(coordinate_type),
      "plant_genetic_georeferenced_source_rows":geocoded,
      "missing_study_ids_among_geocoded":missingids,
      "unresolved_source_binomials_among_geocoded":invalidspecies,
      "distinct_plant_binomials":len(plants),
      "distinct_original_genus_count":len(genera),
      "distinct_original_source_studies":len(original_studies),
      "distinct_original_system_labels":len(systems),
      "distinct_source_pair_geographies":len(source_distinct_coord),
      "repeated_study_species_keys":sum(n>1 for n in group_duplicate.values()),
      "design_metadata_counts":safe_counts(design),
      "experimental_metadata_counts":safe_counts(experimental),
      "experiment_type_metadata_counts":safe_counts(exp_type),
      "trait_type_metadata_counts":safe_counts(trait),
      "elapsed_time_metadata_presence":dict(time_label),
      "species_and_original_study_joint_holdout":split,
      "source_Genetic_label_does_not_prove_adaptive_fitness":True,
      "current_or_ancestral_habitat_continuity_not_evaluated":True,
      "phenotype_or_genetic_response_rows_values_read_for_analysis":False,
      "historical_climate_extracted":False,
      "actual_effect_size_prediction_performed":False,
      "original_TTF_closed_genetics_unmodified":True,
      "raw_sites_or_taxa_rows_exported":False
    }


def main():
    p=argparse.ArgumentParser();p.add_argument("--output",type=Path,required=True)
    args=p.parse_args()
    c=json.loads(CONTRACT.read_text())
    if c["status"]!="RESPONSE_BLIND_PLANT_SOURCE_ELIGIBILITY_FROZEN_BEFORE_TRAIT_RATE_ANALYSIS":
        raise ValueError("Source-only plant protocol not frozen")
    expected=json.loads(NATIVE_LOCK.read_text())["native_original_transport"]
    out={"schema":"ttf_proceed_v62_original_plant_genetic_change_source_eligibility_v01",
         "phenotype_or_genetic_response_rows_values_read_for_analysis":False,
         "actual_effect_size_prediction_performed":False}
    try:
        src=json.loads(SOURCE_CONTRACT.read_text())
        catalog=validated_get_catalog("https://borealisdata.ca/api/datasets/:persistentId/?persistentId=doi%3A10.5683%2FSP3%2FNXSL3Q")
        verify_frozen_catalog(catalog,src["current_source"])
        b,md5,sha256=fetch_original_datafile_bytes(1024206)
        if len(b)!=expected["observed_size_bytes"] or md5!=expected["observed_md5"] or sha256!=expected["observed_sha256"]:
            raise ValueError("Native original source cannot be byte-verified")
        out.update(source_only_plant_counts(b,c))
        out["source_original_md5"]=md5;out["source_original_sha256"]=sha256
    except Exception as e:
        out.update({"status":"PLANT_SOURCE_QUALIFICATION_STOP_NOT_BIOLOGICAL_NULL",
                    "error_type":type(e).__name__,"error_message":str(e)[:350]})
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(out,sort_keys=True,indent=2)+"\n")
    print(json.dumps({k:out.get(k) for k in ("status","plant_genetic_georeferenced_source_rows",
         "distinct_plant_binomials","distinct_original_source_studies",
         "species_and_original_study_joint_holdout")},sort_keys=True))


if __name__=="__main__":main()
