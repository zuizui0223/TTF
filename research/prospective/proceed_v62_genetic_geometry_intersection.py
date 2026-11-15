"""PROCEED v6.2 source-only cross-tabulation of genphen / animal / coordinates.

Full native original CSV is verified by MD5, SHA256 and dataset catalog.
ONLY declared source metadata cells are used. No rate, trait mean, genetic
effect magnitude, historical climate predictor, p-value or AUC is computed.
"""
from __future__ import annotations
import argparse,csv,hashlib,io,json
from collections import Counter,defaultdict
from pathlib import Path

from proceed_v62_borealis_header_only import validated_get_catalog,verify_frozen_catalog
from proceed_v62_metadata_row_census import (
    fetch_original_datafile_bytes,valid_binomial,valid_coordinates,
    safe_counts,connected_components
)

HERE=Path(__file__).resolve().parent
CONTRACT=HERE/"proceed_v62_genetic_geometry_intersection_contract_v01.json"
ROW_CONTRACT=HERE/"proceed_v62_metadata_row_census_contract_v01.json"
SOURCE_CONTRACT=HERE/"proceed_v62_source_header_contract_v01.json"
NATIVE_LOCK=HERE/"proceed_v62_native_original_size_correction_v02.json"


def summarise_stratum(rows):
    taxa={r["species"] for r in rows if r["species"]}
    studies={r["study"] for r in rows if r["study"]}
    genera={s.split()[0] for s in taxa}
    roles={"train":0,"heldout":0}
    graph=defaultdict(set)
    for s in taxa:
        d=hashlib.sha256(("ttf-proceed-v62-genetic-geometry-v01|"+s.lower()).encode()).digest()
        roles["train" if d[0]&1 else "heldout"]+=1
    for r in rows:
        if r["species"] and r["study"]:
            graph[r["study"]].add(r["species"])
    crossing=0
    for group in graph.values():
        role_flags=set()
        for s in group:
            d=hashlib.sha256(("ttf-proceed-v62-genetic-geometry-v01|"+s.lower()).encode()).digest()
            role_flags.add("train" if d[0]&1 else "heldout")
        if len(role_flags)==2:crossing+=1
    components=connected_components(graph)
    return {
        "source_effect_rows":len(rows),
        "distinct_original_species_binomials":len(taxa),
        "distinct_source_studies":len(studies),
        "distinct_genera":len(genera),
        "species_hash_split":roles,
        "studies_crossing_species_hash_roles":crossing,
        "study_species_connected_component_count":len(components),
        "largest_study_species_component":components[0] if components else {"species_count":0,"study_count":0},
    }


def intersection_metadata_summary(raw,design):
    allowed=json.loads(ROW_CONTRACT.read_text())["source_metadata_columns_authorized"]
    if set(allowed)&set(json.loads(ROW_CONTRACT.read_text())["outcome_fields_must_not_be_decoded"]):
        raise ValueError("Outcome columns leaked into data qualification")
    reader=csv.reader(io.TextIOWrapper(io.BytesIO(raw),encoding="utf-8-sig",newline=""),
                      delimiter=",",strict=True)
    header=next(reader)
    if len(header)!=69 or len(set(header))!=69:
        raise ValueError("Source header no longer matches exact 69-column native CSV")
    indices={n:header.index(n) for n in allowed}
    groups=defaultdict(Counter)
    genetic_animal_geo=[];phenotypic_animal_geo=[]
    all_geo_by_genphen=defaultdict(list)
    taxa=Counter();taxa_group_genetic=Counter()
    source_rows=0
    for values in reader:
        if not values:continue
        source_rows+=1
        if len(values)!=69:raise ValueError("Unexpected row-width in exact original source")
        m={k:values[ix].strip() for k,ix in indices.items()}
        original_class=m["genphen"]
        if original_class not in design["genphen_labels_to_count"]:
            raise ValueError("Unexpected source genphen category before outcome opening")
        p1=valid_coordinates(m["sample1_latitude"],m["sample1_longitude"])
        p2=valid_coordinates(m["sample2_latitude"],m["sample2_longitude"])
        geo=("both" if p1 and p2 else
             "first_only" if p1 else "second_only" if p2 else "neither")
        phylum=m["phylum"]
        animal=phylum in design["original_animal_phyla"]
        groups["genphen_by_geostatus"][original_class+"|"+geo]+=1
        groups["phylum_by_genphen_geo"][phylum+"|"+original_class+"|"+geo]+=1
        if p1 and p2:
            species=valid_binomial(m["sp_pub"]) or valid_binomial(m["sp_ncbi"])
            row={
                "species":species,"study":m["ref_ID"],
                "genphen":original_class,"design":m["design"],
                "experimental":m["experimental"],"taxa":m["taxa"],
                "phylum":phylum,
            }
            all_geo_by_genphen[original_class].append(row)
            if animal and original_class=="Genetic":
                genetic_animal_geo.append(row)
                taxa_group_genetic[m["taxa"]]+=1
                groups["genetic_animal_geo_design"][m["design"]]+=1
                groups["genetic_animal_geo_experimental"][m["experimental"]]+=1
                groups["genetic_animal_geo_exp_type"][m["exp_type"]]+=1
                groups["genetic_animal_geo_trait_type"][m["trait_type"]]+=1
                groups["genetic_animal_geo_phylum"][phylum]+=1
            elif animal and original_class=="Phenotypic":
                phenotypic_animal_geo.append(row)
        taxa[m["taxa"]]+=1
    genetic_summary=summarise_stratum(genetic_animal_geo)
    n_species=genetic_summary["distinct_original_species_binomials"]
    n_studies=genetic_summary["distinct_source_studies"]
    genera=genetic_summary["distinct_genera"]
    arms=genetic_summary["species_hash_split"]
    q=design["qualification"]
    if (n_species>=q["min_animal_georeferenced_genetic_label_species"]
        and n_studies>=q["min_source_studies"] and genera>=q["min_genera"]
        and min(arms.values())>=q["min_species_per_future_train_test_arm"]):
        decision="ANIMAL_GENETIC_GEO_DESIGN_SUPPORT_POSSIBLE_HABITAT_AND_HERITABILITY_CONFIRMATION_REQUIRED"
    else:
        decision="ANIMAL_GENETIC_GEO_DESIGN_INSUFFICIENT_NO_RESPONSE_OPENING"
    return {
        "schema":"ttf_proceed_v62_metadata_only_genetic_phenotypic_geolocated_intersection_v01",
        "status":decision,
        "source_rows":source_rows,
        "genphen_by_geostatus":safe_counts(groups["genphen_by_geostatus"]),
        "phylum_by_genphen_geostatus":safe_counts(groups["phylum_by_genphen_geo"],max_labels=120),
        "all_geolocated_genetic_records":summarise_stratum(all_geo_by_genphen["Genetic"]),
        "all_geolocated_phenotypic_records":summarise_stratum(all_geo_by_genphen["Phenotypic"]),
        "animal_phylum_geolocated_genetic":genetic_summary,
        "animal_phylum_geolocated_phenotypic":summarise_stratum(phenotypic_animal_geo),
        "animal_genetic_geolocated_source_design_counts":{
            k:safe_counts(v,max_labels=60) for k,v in groups.items()
            if k.startswith("genetic_animal_geo_")},
        "taxa_label_counts_all":safe_counts(taxa,max_labels=60),
        "taxa_label_counts_animal_genetic_geolocated":safe_counts(taxa_group_genetic,max_labels=60),
        "animal_phylum_is_not_verified_terrestrial_habitat":True,
        "genphen_Genetic_is_source_design_not_adaptive_fitness_proof":True,
        "study_species_cross_arm_leakage_diagnosed":True,
        "no_rate_or_genetic_effect_values_used":True,
        "historical_climate_extracted":False,
        "no_prior_TTF_response_opened":True,
        "public_raw_coordinates_or_rate_data_exported":False
    }


def run():
    p=argparse.ArgumentParser();p.add_argument("--output",required=True,type=Path)
    args=p.parse_args()
    d=json.loads(CONTRACT.read_text())
    if d.get("status")!="SOURCE_RESPONSE_BLIND_STRATIFICATION_FROZEN_BEFORE_RATE_OR_EFFECT_ACCESS":
        raise ValueError("No source-only stratification authority")
    o=json.loads(NATIVE_LOCK.read_text())["native_original_transport"]
    c=json.loads(SOURCE_CONTRACT.read_text())
    out={"schema":"ttf_proceed_v62_metadata_only_genetic_phenotypic_geolocated_intersection_v01",
         "no_rate_or_genetic_effect_values_used":True}
    try:
        catalog=validated_get_catalog("https://borealisdata.ca/api/datasets/:persistentId/?persistentId=doi%3A10.5683%2FSP3%2FNXSL3Q")
        verify_frozen_catalog(catalog,c["current_source"])
        raw,md5,sha256=fetch_original_datafile_bytes(1024206)
        if len(raw)!=o["observed_size_bytes"] or md5!=o["observed_md5"] or sha256!=o["observed_sha256"]:
            raise ValueError("Original native source byte identity changed; stop")
        out.update(intersection_metadata_summary(raw,d))
        out.update({"source_size_bytes":len(raw),"source_md5":md5,"source_sha256":sha256,
                    "response_biological_effects_opened":False})
    except Exception as exc:
        out.update({"status":"SOURCE_META_INTERSECTION_STOP_NOT_BIOLOGICAL_NULL",
                    "error_type":type(exc).__name__,
                    "error_message":str(exc)[:350],
                    "response_biological_effects_opened":False})
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(out,sort_keys=True,indent=2)+"\n")
    print(json.dumps({k:out.get(k) for k in ("status","source_rows",
          "animal_phylum_geolocated_genetic","animal_phylum_geolocated_phenotypic")},
          sort_keys=True))


if __name__=="__main__":run()
