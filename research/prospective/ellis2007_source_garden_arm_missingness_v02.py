"""Separate source-only Ellis original experimental-arm missingness audit.

Original 2007 primary field adult transplant: each of six gardens received
four sources, one local, one other quartz habitat within its basin, and two
same-quartz sources in other basins. Hargreaves secondary original has only
three unique source IDs per garden. Count missing classes from exact original
label strings, without touching numerical fitness/phenotype/genetic outcomes.

A missing Hargreaves source×garden edge means missing from this curated
secondary extract, NOT absence from the primary experimental cohort.
"""
from __future__ import annotations
import argparse
import hashlib
import json
from collections import Counter,defaultdict
from pathlib import Path
from bontrager_hargreaves_finegrain_overlap_v01 import original_bytes,haversine
from hargreaves_source_as_garden_georecovery_v01 import source_metadata,one_point
from ellis2007_six_sites_basin_habitat_geometry_v01 import source_id_category

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
CONTRACT=HERE/"ellis2007_missing_transplant_arms_diagnosis_contract_v02.json"

CATEGORIES=(
    "local_same_source_and_garden",
    "same_basin_other_quartz",
    "other_basin_same_quartz",
    "other_basin_other_quartz_not_in_paper_design",
)

def classify(source:dict,garden:dict) -> str:
    if not (source.get("basin") and source.get("habitat")
            and garden.get("basin") and garden.get("habitat")):
        raise ValueError("Original site source/garden label did not encode unique basin and quartz type")
    same_basin=source["basin"]==garden["basin"]
    same_quartz=source["habitat"]==garden["habitat"]
    return (
        CATEGORIES[0] if same_basin and same_quartz else
        CATEGORIES[1] if same_basin else
        CATEGORIES[2] if same_quartz else CATEGORIES[3]
    )

def arm_census(meta:dict,contract:dict)->dict:
    if contract["status"]!="FROZEN_AFTER_ORIGINAL_SOURCE_ONLY_18_OF_24_DESIGN_DISCOVERY_BEFORE_CATEGORY_BREAKDOWN":
        raise ValueError("Ellis arm classifier not precommitted")
    primary=contract["predeclared_primary_article_transplant_topology"]
    if (primary["expected_total_pair_edges"]!=24
        or primary["expected_local_same_site_pairs"]!=6
        or primary["expected_same_basin_different_quartz_pairs"]!=6
        or primary["expected_other_basin_same_quartz_pairs"]!=12
        or primary["expected_other_basin_different_quartz_pairs"]!=0):
        raise ValueError("Source-only primary original expected topology changed")
    original=contract["previous_sha_verified_geometry"]
    study=original["study_label"];taxon=original["species"]
    origin=sorted(k for k in meta["origin_keys"] if k[0]==study and k[1]==taxon)
    if len(origin)!=6:raise ValueError("Original source population count drift")
    code={k[2]:source_id_category(k[2],{
        "study":{"published_original_methods":{
          "focal_drainage_basins":["Ariz","Sout","QK"],
          "quartz_habitat_groups":["dense","sparse"]
        }}
    }) for k in origin}
    if any(not c["basin"] or not c["habitat"] for c in code.values()):
        raise ValueError("Not all original source labels are directly basin×quartz classifiable")
    if len({(v["basin"],v["habitat"]) for v in code.values()})!=6:
        raise ValueError("Source labels not a unique 3×2 habitat-basin Cartesian grid")
    sites={k[2]:one_point(meta["gardens"].get(k,set())) for k in origin}
    if any(v is None for v in sites.values()) or len(set(sites.values()))!=6:
        raise ValueError("Original garden-derived locations ambiguous or missing")

    per_garden=defaultdict(set)
    # Earlier study group join counts each exact original source-garden pair
    # only once across experiment labels and fitness-TYPE labels; never
    # model original survival/growth/reproduction numbers.
    for cell,sources in meta["source_groups"].items():
        if cell[0]!=study or cell[1]!=taxon:continue
        garden=cell[3]
        if garden not in code:
            raise ValueError("Hargreaves observation uses unrecognized original garden ID")
        for source in sources:
            if source not in code:
                raise ValueError("Hargreaves observation uses unrecognized original source ID")
            per_garden[garden].add(source)
    if len(per_garden)!=6:raise ValueError("Original source-only garden count drift")
    observed_edges={(source,garden) for garden,sources in per_garden.items() for source in sources}
    if (len(observed_edges)!=original["observed_origin_destination_pairs"]
        or set(len(v) for v in per_garden.values())!={3}):
        raise ValueError("Prior first-run original source-garden arm geometry drift")

    expected=set()
    for garden,site_info in code.items():
        for source,from_info in code.items():
            if classify(from_info,site_info)!=CATEGORIES[3]:
                expected.add((source,garden))
    if len(expected)!=24 or not observed_edges.issubset(expected):
        raise ValueError("Observed secondary source edges not compatible with published design")
    missing=expected-observed_edges
    count_obs=Counter(classify(code[src],code[dst]) for src,dst in observed_edges)
    count_missing=Counter(classify(code[src],code[dst]) for src,dst in missing)
    count_expected=Counter(classify(code[src],code[dst]) for src,dst in expected)
    orig={
        "local_same_source_and_garden":primary["expected_local_same_site_pairs"],
        "same_basin_other_quartz":primary["expected_same_basin_different_quartz_pairs"],
        "other_basin_same_quartz":primary["expected_other_basin_same_quartz_pairs"],
        "other_basin_other_quartz_not_in_paper_design":primary["expected_other_basin_different_quartz_pairs"]
    }
    if any(count_expected[cat]!=orig[cat] for cat in CATEGORIES):
        raise AssertionError("Original primary design classification inconsistent with 3×2 grid")
    basin_dists={}
    for basin in ("Ariz","Sout","QK"):
        k=[name for name in sites if code[name]["basin"]==basin]
        if len(k)!=2:raise AssertionError("Exactly two original sites per basin required")
        basin_dists[basin]=round(haversine(sites[k[0]],sites[k[1]]),5)

    return {
      "schema":"ttf_Ellis2007_primary_24_source_garden_vs_Hargreaves_18_source_only_edges_diagnosis_v02",
      "status":"PRIMARY_ORIGINAL_24_SOURCE_GARDEN_DESIGN_VS_SECONDARY_EXTRACT_18_WITH_SOURCE_ARM_MISSINGNESS",
      "original_paper_doi":primary["DOI"],
      "Hargreaves_original_study_label":study,
      "Hargreaves_original_SHA256_verified":original["source_original_hargreaves_sha256"],
      "same_2000_original_adult_transplant_six_populations":6,
      "original_primary_expected_directed_source_garden_pairs":len(expected),
      "observed_Hargreaves_distinct_directed_source_garden_pairs":len(observed_edges),
      "source_garden_pairs_absent_from_Hargreaves_extract":len(missing),
      "published_origin_site_pairs_by_class":[
          {"class":cat,"paper_expected":count_expected[cat],
           "original_Hargreaves_extract_observed":count_obs[cat],
           "not_represented_in_secondary_original":count_missing[cat]}
          for cat in CATEGORIES
      ],
      "published_source_population_IDs_per_target_garden":4,
      "observed_secondary_source_IDs_per_target_garden":3,
      "original_basin_label_same_quartz_source_pair_verified_from_ID":True,
      "within_basin_dense_sparse_pair_distance_km":basin_dists,
      "publisher_within_basin_approximate_distance_claim_km":[0.4,1.2],
      "within_basin_pairs_beyond_1p2km":sum(v>1.2 for v in basin_dists.values()),
      "coordinate_discrepancy_cause_certified":False,
      "source_pair_absence_is_not_original_2000_intervention_absence":True,
      "source_pair_absence_may_be_downstream_meta_analysis_extraction":True,
      "raw_original_location_or_plant_ID_values_exported":False,
      "no_fitness_growth_survival_reproduction_or_AFLP_response_rows_read":True,
      "no_LGM_21ka_0BP_climate_effect_estimated":True,
      "old_Bontrager_44of50_terminal_stop_unchanged":True
    }

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--output",required=True,type=Path)
    args=p.parse_args()
    c=json.loads(CONTRACT.read_text(encoding="utf-8"))
    src=json.loads((ROOT/c["source_metadata"]["original_hargreaves_archive_contract"]).read_text())
    raw=original_bytes(src["original"]["hargreaves"])
    if hashlib.sha256(raw).hexdigest()!=c["previous_sha_verified_geometry"]["source_original_hargreaves_sha256"]:
        raise ValueError("Original 2007 secondary source original SHA changed")
    rule=json.loads((ROOT/c["source_metadata"]["read_metadata_helper_contract"]).read_text())
    meta=source_metadata(raw,"hargreaves",rule)
    x=arm_census(meta,c)
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(x,sort_keys=True,indent=2,ensure_ascii=False)+"\n")
    print("ELLIS2007_SOURCE_EXTRACT_MISSING_ARM_CATEGORY=",json.dumps({
      "original_paper_expected":x["original_primary_expected_directed_source_garden_pairs"],
      "secondary_observed":x["observed_Hargreaves_distinct_directed_source_garden_pairs"],
      "missing":x["source_garden_pairs_absent_from_Hargreaves_extract"],
      "by_original_design_arm":x["published_origin_site_pairs_by_class"],
      "within_basin_three_distances":x["within_basin_dense_sparse_pair_distance_km"],
      "fitness_outcome_unopened":x["no_fitness_growth_survival_reproduction_or_AFLP_response_rows_read"]
    },sort_keys=True))

if __name__=="__main__":
    main()
