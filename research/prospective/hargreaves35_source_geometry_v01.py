"""Hargreaves's 35 unlinked original-study labels: source-only feasibility audit.

Read from official SHA/MD5-pinned Zenodo original CSV. The original CSV lacks
origin population lat/lon; source IDs and garden IDs do not authorize
historical exposure or experimental independence inference. No fitness,
treatment effect, viability, trait response or effect-size field is indexed.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
from collections import Counter, defaultdict
from pathlib import Path

from bontrager_hargreaves_finegrain_overlap_v01 import (
    norm_taxon, valid_coord, original_bytes,
)

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
CONTRACT=HERE/"hargreaves35_source_geometry_contract_v01.json"
PLANT_GROUPS=("annual","herbaceous_perennial","woody_perennial")
ALLOWED=("study","functional_group","taxon","source","site","experiment",
         "fitness_type","year_begin","lat_site","long_site")


def audit_hargreaves35(raw:bytes, contract:dict, receipt:dict) -> dict:
    src=contract["source"]
    if (len(raw)!=src["exact_original_byte_count"]
        or hashlib.md5(raw).hexdigest()!=src["expected_original_md5"]
        or hashlib.sha256(raw).hexdigest()!=src["expected_original_sha256"]):
        raise ValueError("Original Hargreaves bytes/digests no longer match prior source")
    if hashlib.sha256(raw.split(b"\n",1)[0]+b"\n").hexdigest()!=src["expected_header_sha256"]:
        raise ValueError("Original first-line checksum drifted")
    if tuple(contract["allowed_columns"])!=ALLOWED or "fitness" in ALLOWED:
        raise ValueError("Response-independent allowlist drift")
    if tuple(contract["group_allowlist"])!=PLANT_GROUPS:
        raise ValueError("Original plant functional group eligibility drift")
    targets=tuple(contract["reference"]["exact_unlinked_study_labels"])
    if len(targets)!=35 or len(set(targets))!=35 or sorted(targets)!=sorted(receipt["Hargreaves_unlinked_full_study_labels_for_followup"]):
        raise ValueError("Prior exact 35 unlinked study labels changed")
    if receipt["Hargreaves_study_labels_not_in_any_Bontrager_original_plant_row"]!=35:
        raise ValueError("Prior first-shot source-only study counts not reproduced")

    if src["expected_text_encoding"]!="cp1252":
        raise ValueError("Published source-encoding correction not authorized")
    txt=raw.decode("cp1252",errors="strict")
    reader=csv.reader(io.StringIO(txt,newline=""),strict=True)
    header=next(reader)
    if len(header)!=47 or len(set(header))!=47 or not set(ALLOWED).issubset(header):
        raise ValueError("Original Hargreaves 47-column schema or required metadata drift")
    if any(field in header for field in src["source_lonlat_columns_absent"]):
        raise ValueError("New source location columns present: require new prospectively frozen contract")
    indices={field:header.index(field) for field in ALLOWED}
    refs={key:set() for key in targets}
    cells=defaultdict(set)
    gardens=defaultdict(set)
    experiment_tokens=defaultdict(set)
    garden_tokens=defaultdict(set)
    source_tokens=defaultdict(set)
    meta_study=defaultdict(Counter)
    strata_year=defaultdict(set)
    st_taxa=defaultdict(set)
    missing=Counter()
    n_original=0
    for row in reader:
        if not row:
            continue
        n_original+=1
        if n_original>200000 or len(row)!=47:
            raise ValueError("Source row width or bounded original row count drift")
        # Source metadata only; never look up "fitness" or response columns.
        study=row[indices["study"]].strip()
        if study not in refs:
            continue
        group=row[indices["functional_group"]].strip()
        if group not in PLANT_GROUPS:
            continue
        meta_study[study]["source_only_plant_rows"]+=1
        taxon=norm_taxon(row[indices["taxon"]])
        source=row[indices["source"]].strip()
        site=row[indices["site"]].strip()
        experiment=row[indices["experiment"]].strip()
        ftype=row[indices["fitness_type"]].strip()
        start_year=row[indices["year_begin"]].strip()
        dest=valid_coord(row[indices["long_site"]],row[indices["lat_site"]])
        if taxon:
            st_taxa[study].add(taxon)
        else:
            missing["nonbinomial_taxon"]+=1
        if source:
            source_tokens[study].add(source)
        else:
            missing["empty_source_ID"]+=1
        if site:
            garden_tokens[study].add(site)
        else:
            missing["empty_garden_ID"]+=1
        if experiment:
            experiment_tokens[study].add(experiment)
        else:
            missing["missing_original_experiment_label"]+=1
        if start_year:
            strata_year[study].add(start_year)
        if not (taxon and source and site and ftype):
            missing["incomplete_source_garden_metric_key"]+=1
            continue
        key=(study,taxon,experiment if experiment else "<UNKNOWN_EXPERIMENT>",site,ftype)
        cells[key].add(source)
        gardens[key].add(dest)
        refs[study].add(key)
    if n_original!=11088:
        raise ValueError("Hargreaves original 11088 physical rows drifted")
    per_study=[]
    for study in sorted(targets):
        groups=refs[study]
        qualifying=[k for k in groups if len(cells[k])>=2]
        count_geo=sum(len(gardens[k])==1 and None not in gardens[k] for k in qualifying)
        count_missing_experiment=sum(k[2]=="<UNKNOWN_EXPERIMENT>" for k in qualifying)
        per_study.append({
            "original_study_label":study,
            "source_only_plant_rows":meta_study[study]["source_only_plant_rows"],
            "candidate_original_plant_taxa":len(st_taxa[study]),
            "distinct_original_source_IDs":len(source_tokens[study]),
            "distinct_original_garden_IDs":len(garden_tokens[study]),
            "distinct_experiment_labels":len(experiment_tokens[study]),
            "distinct_start_year_labels":len(strata_year[study]),
            "candidate_comparable_same_garden_cells":len(qualifying),
            "candidate_same_garden_cells_with_complete_garden_coordinate":count_geo,
            "candidate_same_garden_cells_missing_or_ambiguous_garden_coordinate":len(qualifying)-count_geo,
            "candidate_same_garden_cells_missing_experiment_label":count_missing_experiment,
            "source_origin_geocoordinates_available_in_original":False,
            "original_primary_study_experiment_independence_certified":False,
            "does_not_requalify_stopped_adaptation_macro_panel":True,
        })
    summary={
        "original_study_labels_previously_unlinked":len(per_study),
        "study_labels_with_at_least_one_plant_row":sum(x["source_only_plant_rows"]>0 for x in per_study),
        "study_labels_with_at_least_two_source_IDs":sum(x["distinct_original_source_IDs"]>=2 for x in per_study),
        "study_labels_with_same_garden_multi_source_opportunity":sum(x["candidate_comparable_same_garden_cells"]>0 for x in per_study),
        "sum_candidate_same_garden_cells_not_independent_experiments":sum(x["candidate_comparable_same_garden_cells"] for x in per_study),
        "sum_candidate_cells_missing_experiment_ID":sum(x["candidate_same_garden_cells_missing_experiment_label"] for x in per_study),
        "study_labels_with_no_original_experiment_ID":sum(x["distinct_experiment_labels"]==0 for x in per_study),
        "study_labels_with_multiple_original_garden_IDs":sum(x["distinct_original_garden_IDs"]>=2 for x in per_study),
        "plant_rows_in_35_unlinked_source_labels":sum(x["source_only_plant_rows"] for x in per_study),
        "verified_source_population_georeferenced_studies_in_Hargreaves_original":0,
        "verified_independent_original_experiments":0,
        "fitness_response_values_read":0,
        "no_new_historical_predictive_fitness_model":True,
    }
    return {
        "schema":"ttf_Hargreaves35_metadata_only_unlinked_source_garden_opportunity_v01",
        "status":"SOURCE_ONLY_CANDIDATE_GEOMETRY_NO_ORIGIN_COORDINATES_NO_INDEPENDENCE_CERTIFICATION",
        "source_doi":src["data_doi"],
        "source_SHA256_verified":src["expected_original_sha256"],
        "frozen_unlinked_source_study_list_receipt_SHA256":contract["reference"]["metadata_receipt_original_sha256"],
        "original_CSV_header_has_source_lat_lon":False,
        "metadata_allowlist":list(ALLOWED),
        "all_source_row_count":n_original,
        "missing_metadata_counts":dict(sorted(missing.items())),
        "summary":summary,
        "studies":per_study,
        "fitness_values_or_effects_indexed":False,
        "Bontrager_44_of_50_frozen_stop_unchanged":True,
        "Hargreaves_35_labels_not_promoted_to_independent_studies":True,
        "raw_source_or_garden_IDs_or_coordinates_exported":False,
    }


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("--output",type=Path,required=True)
    args=parser.parse_args()
    c=json.loads(CONTRACT.read_text(encoding="utf-8"))
    if c.get("status")!="FROZEN_BEFORE_HARGREAVES35_METADATA_ROW_CENSUS":
        raise ValueError("No source-only predata authority for Hargreaves35")
    receipt_path=ROOT / c["reference"]["metadata_receipt"]
    receipt_bytes=receipt_path.read_bytes()
    if hashlib.sha256(receipt_bytes).hexdigest()!=c["reference"]["metadata_receipt_original_sha256"]:
        raise ValueError("Frozen 35-study prior source-only receipt bytes drifted")
    receipt=json.loads(receipt_bytes)
    source=c["source"]
    entry={"doi":source["data_doi"],"zenodo_record":source["zenodo_record_id"],
           "file":source["file_name"],"size":source["exact_original_byte_count"],
           "md5":source["expected_original_md5"],"sha256":source["expected_original_sha256"]}
    original=original_bytes(entry)
    out=audit_hargreaves35(original,c,receipt)
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(out,ensure_ascii=False,sort_keys=True,indent=2)+"\n",encoding="utf-8")
    print("HARGREAVES35_SOURCE_GEOMETRY=",json.dumps(out["summary"],sort_keys=True))


if __name__=="__main__":
    main()
