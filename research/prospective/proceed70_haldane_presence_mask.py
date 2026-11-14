"""PROCEED70 response-size-only *presence flags*, no Haldane arithmetic or rates.

The source study/species/site universe is frozen on response-BLIND geography.
For fields mean1/mean2/SD/n, this stage inspects only whether their raw CSV
text is absent, never converts them to numbers or evaluates magnitudes.
Read code provenance and original source MD5/SHA. No GENETIC/phenotypic
association or any previously frozen TTF outcome is opened.
"""
from __future__ import annotations
import argparse,csv,io,json
from collections import Counter,defaultdict
from pathlib import Path
from proceed_v62_borealis_header_only import validated_get_catalog,verify_frozen_catalog
from proceed_v62_metadata_row_census import (fetch_original_datafile_bytes,
                                              valid_binomial,valid_coordinates)
from proceed_plant_chelsa_source_support import source_geometry,component_assignment
from proceed_v62_genetic_geometry_intersection import SOURCE_CONTRACT,NATIVE_LOCK

HERE=Path(__file__).resolve().parent
CONTRACT=HERE/"proceed70_haldane_presence_mask_contract_v01.json"
FIELDS=("mean1","mean2","sd1","sd2","n1","n2","generations")
META=("ref_ID","sp_pub","sp_ncbi","phylum","taxa","genphen",
      "sample1_latitude","sample1_longitude",
      "sample2_latitude","sample2_longitude",
      "trait_type","design")


def is_present_raw_token(value,empty_tokens):
    """Only binary CSV-cell existence: do not parse or compare a magnitude."""
    return str(value).strip().lower() not in empty_tokens


def missingness_summary(raw,protocol):
    fields=protocol["field_presence_only"]["columns"]
    if fields!=list(FIELDS):
        raise ValueError("Frozen Haldane presence-field set changed")
    null={str(z).lower() for z in protocol["field_presence_only"]["absent_if_trimmed_lowercase_in"]}
    if any(z in META for z in FIELDS):
        raise ValueError("Response presence fields mixed with locality metadata")
    metadata_reader=csv.reader(io.TextIOWrapper(io.BytesIO(raw),encoding="utf-8-sig",newline=""),
                              delimiter=",",strict=True)
    header=next(metadata_reader,[])
    if len(header)!=69 or len(set(header))!=len(header):
        raise ValueError("Unrecognized PROCEED v6.2 native 69-field CSV")
    if not set(META+FIELDS).issubset(header):
        raise ValueError("Source metadata or required presence field absent")
    mi={k:header.index(k) for k in META}
    ri={k:header.index(k) for k in FIELDS}
    rows=0;all_source=0;qual_rows=0
    source_study_species=defaultdict(set)
    encountered_mask=Counter();present_fields=Counter()
    usable_geom=set();overall_geom=set();usable_studies=set();usable_species=set()
    usability_by_trait=Counter();usability_by_design=Counter()
    for values in metadata_reader:
        if not values:continue
        all_source+=1
        if len(values)!=69:raise ValueError("Original source row width changed")
        m={name:values[ix].strip() for name,ix in mi.items()}
        if m["phylum"]!="Streptophyta" or m["taxa"]!="Plant" or m["genphen"]!="Genetic":
            continue
        c1=valid_coordinates(m["sample1_latitude"],m["sample1_longitude"])
        c2=valid_coordinates(m["sample2_latitude"],m["sample2_longitude"])
        taxon=valid_binomial(m["sp_pub"]) or valid_binomial(m["sp_ncbi"])
        study=m["ref_ID"]
        if not(c1 and c2 and taxon and study):
            continue
        rows+=1
        source_study_species[study].add(taxon)
        geom=(study,taxon,c1,c2)
        overall_geom.add(geom)
        # We NEVER decode or convert a rate, mean, sd or sample size value:
        # only a per-column missingness bit for a *predeclared* set.
        indicators=tuple(is_present_raw_token(values[ri[name]],null) for name in FIELDS)
        for name,ok in zip(FIELDS,indicators):
            if ok:present_fields[name]+=1
        encountered_mask["".join("1" if z else "0" for z in indicators)]+=1
        if all(indicators):
            qual_rows+=1;usable_geom.add(geom)
            usable_studies.add(study);usable_species.add(taxon)
            usability_by_trait[m["trait_type"] or "<missing>"]+=1
            usability_by_design[m["design"] or "<missing>"]+=1

    scope=protocol["existing_scope"]
    if (all_source!=9263 or rows!=scope["original_source_rows"] or
        len(source_study_species)!=scope["original_studies"] or
        len(overall_geom)!=scope["original_pair_geometries"]):
        raise ValueError("Original 70-plant source-only geometry changed")
    roles,split,digest=component_assignment(source_study_species)
    if (digest!=scope["frozen_study_species_split_sha256"] or
        split["species_per_role"]!={"train":35,"heldout":35}):
        raise ValueError("Frozen 35/35 original study-species separation changed")
    mask_species_roles={"train":set(),"heldout":set()}
    geom_role_counts=Counter()
    for study,taxon,c1,c2 in usable_geom:
        role=roles[taxon]
        mask_species_roles[role].add(taxon)
        geom_role_counts[role]+=1
    groups={x.split()[0] for x in usable_species}
    bounds=protocol["scope_support_criteria"]
    enough=(len(usable_species)>=bounds["min_species_total"]
            and len(usable_studies)>=bounds["min_studies_total"]
            and len(groups)>=bounds["min_genera_total"]
            and len(mask_species_roles["train"])>=bounds["min_train_species"]
            and len(mask_species_roles["heldout"])>=bounds["min_heldout_species"]
            and geom_role_counts["train"]>=bounds["min_train_distinct_study_species_pair_geometries"]
            and geom_role_counts["heldout"]>=bounds["min_heldout_distinct_study_species_pair_geometries"])
    return {
        "schema":"ttf_PROCEED70_source_only_binary_Haldane_component_presence_v01",
        "status":bounds["pass_status"] if enough else bounds["fail_status"],
        "all_source_rows":all_source,
        "original_geographic_plant_source_rows":rows,
        "original_geographic_study_species_pairs":len(overall_geom),
        "rows_with_all_required_Haldane_source_fields_present":qual_rows,
        "distinct_original_study_species_geographic_pairs_with_any_complete_mask":len(usable_geom),
        "distinct_species_with_any_complete_mask":len(usable_species),
        "distinct_studies_with_any_complete_mask":len(usable_studies),
        "distinct_genera_with_any_complete_mask":len(groups),
        "complete_mask_original_species_by_role":{k:len(v) for k,v in mask_species_roles.items()},
        "complete_mask_unique_pair_geometries_by_role":dict(sorted(geom_role_counts.items())),
        "source_presence_only_counts_per_required_field":dict(sorted(present_fields.items())),
        "source_presence_bitmask_counts_frozen_field_order":dict(sorted(encountered_mask.items())),
        "source_trait_type_counts_on_complete_mask":dict(sorted(usability_by_trait.items())),
        "source_study_design_counts_on_complete_mask":dict(sorted(usability_by_design.items())),
        "original_study_species_role_sha256":digest,
        "numeric_mean_SD_n_rate_or_fitness_values_decoded":False,
        "Haldane_or_Darwin_rates_calculated":False,
        "source_numeric_eligibility_and_outlier_QC_not_performed":True,
        "genetic_trait_adaptation_or_plasticity_prediction_not_tested":True,
        "previous_TTF_genetic_outcomes_unchanged":True,
        "original_sites_and_phenotype_values_not_exported":True,
    }


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--output",required=True,type=Path)
    args=p.parse_args()
    contract=json.loads(CONTRACT.read_text())
    if contract["status"]!="FROZEN_BEFORE_ANY_NUMERIC_TRAIT_OR_RATE_VALUE_ACCESS":
        raise ValueError("No frozen original Haldane presence-mask authority")
    src=json.loads(SOURCE_CONTRACT.read_text())["current_source"]
    native=json.loads(NATIVE_LOCK.read_text())["native_original_transport"]
    output={"schema":"ttf_PROCEED70_source_only_binary_Haldane_component_presence_v01",
            "numeric_mean_SD_n_rate_or_fitness_values_decoded":False,
            "Haldane_or_Darwin_rates_calculated":False,
            "genetic_trait_adaptation_or_plasticity_prediction_not_tested":True,
            "previous_TTF_genetic_outcomes_unchanged":True}
    try:
        catalog=validated_get_catalog(
            "https://borealisdata.ca/api/datasets/:persistentId/?persistentId=doi%3A10.5683%2FSP3%2FNXSL3Q"
        )
        verify_frozen_catalog(catalog,src)
        raw,md5,sha=fetch_original_datafile_bytes(1024206)
        if len(raw)!=native["observed_size_bytes"] or md5!=native["observed_md5"] or sha!=native["observed_sha256"]:
            raise ValueError("Original source size MD5 or SHA256 changed before response presence-mask check")
        output.update(missingness_summary(raw,contract))
        output["original_rate_source_sha256"]=sha
        output["original_rate_source_md5"]=md5
    except Exception as exc:
        output["status"]="SOURCE_RESPONSE_PRESENCE_MASK_STOP_NOT_BIOLOGICAL_NULL"
        output["error_type"]=type(exc).__name__
        output["error_message"]=str(exc)[:400]
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(output,sort_keys=True,indent=2)+"\n")
    print(json.dumps({k:output.get(k) for k in (
        "status","rows_with_all_required_Haldane_source_fields_present",
        "distinct_species_with_any_complete_mask",
        "distinct_studies_with_any_complete_mask",
        "complete_mask_original_species_by_role")},sort_keys=True))


if __name__=="__main__":
    main()
