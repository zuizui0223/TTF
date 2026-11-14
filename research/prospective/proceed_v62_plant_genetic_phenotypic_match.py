"""Find *source metadata matches*, not genetic/plasticity effect differences.

Only category labels and original collection site coordinates are used. No
PROCEED numerical rate, standard deviation, genetic value, or fitness result
is accessed, fitted or exported.
"""
from __future__ import annotations
import argparse,csv,io,json
from collections import Counter,defaultdict
from pathlib import Path

from proceed_v62_metadata_row_census import (
    fetch_original_datafile_bytes,valid_binomial,valid_coordinates)
from proceed_v62_borealis_header_only import validated_get_catalog,verify_frozen_catalog
from proceed_v62_genetic_geometry_intersection import SOURCE_CONTRACT,NATIVE_LOCK

HERE=Path(__file__).resolve().parent
CONTRACT=HERE/"proceed_v62_plant_genetic_phenotypic_match_contract_v01.json"


def matched_summary(original, contract):
    allowed=contract["source_metadata_fields_allowed"]
    if set(allowed)&set(contract["forbidden_rate_columns"]):
        raise ValueError("Rate response field in source-only allowlist")
    csv_read=csv.reader(io.TextIOWrapper(io.BytesIO(original),encoding="utf-8-sig",newline=""),
                        delimiter=",",strict=True)
    names=next(csv_read)
    if len(names)!=69 or len(set(names))!=69 or not set(allowed).issubset(names):
        raise ValueError("Original native PROCEED v6.2 file schema changed")
    positions={k:names.index(k) for k in allowed}
    categories=("Genetic","Phenotypic")
    strict=defaultdict(set)
    intermediate=defaultdict(set)
    species_trait=defaultdict(set)
    source_rows=0;plants=0;complete_plants=0;valid_taxa=set()
    label_counts=Counter()
    for vals in csv_read:
        if not vals:continue
        source_rows+=1
        if len(vals)!=69:raise ValueError("Unexpected source row width")
        # Only source study, category, original population coordinates and
        # general trait family are decoded. Original rate values are excluded.
        m={k:vals[i].strip() for k,i in positions.items()}
        if m["phylum"]!="Streptophyta" or m["taxa"]!="Plant":
            continue
        if m["genphen"] not in categories:
            raise ValueError("Unexpected PROCEED source genphen label")
        plants+=1
        a=valid_coordinates(m["sample1_latitude"],m["sample1_longitude"])
        b=valid_coordinates(m["sample2_latitude"],m["sample2_longitude"])
        sp=valid_binomial(m["sp_pub"]) or valid_binomial(m["sp_ncbi"])
        study=m["ref_ID"]
        trait=m["trait_type"]
        if not(a and b and sp and study and trait):
            continue
        complete_plants+=1
        valid_taxa.add(sp)
        label_counts[m["genphen"]]+=1
        strict[(study,sp,trait,a,b)].add(m["genphen"])
        intermediate[(study,sp,trait)].add(m["genphen"])
        species_trait[(sp,trait)].add(m["genphen"])

    good_strict=[key for key,v in strict.items() if len(v)==2]
    good_study=[key for key,v in intermediate.items() if len(v)==2]
    good_species=[key for key,v in species_trait.items() if len(v)==2]
    strict_species=set(v[1] for v in good_strict)
    strict_studies=set(v[0] for v in good_strict)
    c=contract["gates"]
    passed=(len(strict_species)>=c["minimum_strictly_paired_distinct_plant_species"] and
            len(strict_studies)>=c["minimum_strictly_paired_original_studies"])
    return {
        "schema":"ttf_proceed_v62_genetic_phenotypic_plant_source_metadata_exact_matches_v01",
        "status":"SOURCE_PAIRED_DESIGN_POSSIBLE_NEEDS_ORIGINAL_GENOTYPE_AND_FITNESS_AUDIT" if passed else c["if_fail"],
        "all_original_source_rows":source_rows,
        "all_source_plant_genetic_and_phenotypic_rows":plants,
        "plant_valid_geography_taxon_study_trait_metadata_rows":complete_plants,
        "georeferenced_plant_source_category_counts":dict(label_counts),
        "georeferenced_plant_distinct_species_across_both_source_categories":len(valid_taxa),
        "strict_same_study_species_trait_ordered_geography_matched_keys":len(good_strict),
        "strictly_matched_distinct_species":len(strict_species),
        "strictly_matched_distinct_original_studies":len(strict_studies),
        "same_study_species_trait_matched_keys_regardless_of_geography":len(good_study),
        "same_study_species_trait_distinct_species":len(set(x[1] for x in good_study)),
        "same_species_trait_matched_keys_regardless_of_study_geography":len(good_species),
        "same_species_trait_distinct_species":len(set(x[0] for x in good_species)),
        "no_numerical_trait_or_genetic_rate_access":True,
        "source_categories_not_plasticity_and_selection_decomposition":True,
        "no_new_climate_response_association_fitted":True,
        "old_TTF_genetic_outcomes_unchanged":True,
        "no_original_site_coordinates_or_effect_values_exported":True,
    }


def run():
    p=argparse.ArgumentParser();p.add_argument("--output",required=True,type=Path)
    args=p.parse_args()
    c=json.loads(CONTRACT.read_text())
    if c["status"]!="FROZEN_PRIOR_TO_SOURCE_ONLY_GEXP_PEXP_MATCHING_NO_RATE_VALUES_OPENED":
        raise ValueError("Source-only paired Genetic/Phenotypic protocol not frozen")
    native=json.loads(NATIVE_LOCK.read_text())["native_original_transport"]
    src=json.loads(SOURCE_CONTRACT.read_text())["current_source"]
    report={
        "schema":"ttf_proceed_v62_genetic_phenotypic_plant_source_metadata_exact_matches_v01",
        "no_numerical_trait_or_genetic_rate_access":True,
        "no_new_climate_response_association_fitted":True,
        "old_TTF_genetic_outcomes_unchanged":True,
        "no_original_site_coordinates_or_effect_values_exported":True
    }
    try:
        catalog=validated_get_catalog("https://borealisdata.ca/api/datasets/:persistentId/?persistentId=doi%3A10.5683%2FSP3%2FNXSL3Q")
        verify_frozen_catalog(catalog,src)
        raw,md5,sha=fetch_original_datafile_bytes(1024206)
        if (len(raw)!=native["observed_size_bytes"] or
            md5!=native["observed_md5"] or
            sha!=native["observed_sha256"]):
            raise ValueError("Original native source MD5, SHA256 or byte size changed")
        report.update(matched_summary(raw,c))
        report["exact_source_md5"]=md5
        report["exact_source_sha256"]=sha
    except Exception as exc:
        report["status"]="NOT_EVALUABLE_SOURCE_METADATA_PAIRING_OR_REPLAY"
        report["error_type"]=type(exc).__name__
        report["error_message"]=str(exc)[:300]
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(report,sort_keys=True,indent=2)+"\n")
    print(json.dumps({k:report.get(k) for k in (
      "status","strictly_matched_distinct_species",
      "strictly_matched_distinct_original_studies",
      "same_study_species_trait_matched_keys_regardless_of_geography")},
      sort_keys=True))


if __name__=="__main__":run()
