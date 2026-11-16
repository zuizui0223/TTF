"""Ellis et al. (2007) original six-population source-only geometry and design.

Published Methods establish three drainage basins (Ariz, Sout, QK), with dense
and sparse quartz habitats in each. Extract only allowed study/taxon/origin/
garden/fitness-*type* labels and destination coordinates from the exact
SHA-verified Hargreaves original. Never read fitness values, history covariates,
genetic outcomes, or original individual raw observations.

Compute all 15 anonymous pair distances and source×garden comparison topology.
Classify basin and quartz group ONLY when labels explicitly encode unique
names, never from site geography or a chosen matching distance.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import re
from collections import Counter, defaultdict
from pathlib import Path

from bontrager_hargreaves_finegrain_overlap_v01 import original_bytes,haversine
from hargreaves_source_as_garden_georecovery_v01 import (
    source_metadata,one_point,percent_at
)

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
CONTRACT=HERE/"ellis2007_six_sites_basin_habitat_geometry_contract_v01.json"

def source_id_category(value:str,c:dict)->dict:
    """Strict observational lexical match; ambiguous/no match gives null."""
    low=value.casefold()
    r=c["study"]["published_original_methods"]
    basins=[b for b in r["focal_drainage_basins"] if b.casefold() in low]
    habitats=[h for h in r["quartz_habitat_groups"] if h.casefold() in low]
    return {
        "basin":basins[0] if len(basins)==1 else None,
        "habitat":habitats[0] if len(habitats)==1 else None
    }

def bucket(distance:float,limits:list[float])->str:
    for x in limits:
        if distance<=x:
            return "le_"+str(x).replace(".","p")+"_km"
    return "gt_"+str(limits[-1]).replace(".","p")+"_km"

def geometry(meta:dict,c:dict)->dict:
    if c["status"]!="FROZEN_BEFORE_ELLIS2007_ORIGINAL_HARGREAVES_GEOMETRY_AUDIT":
        raise ValueError("Original source-only Ellis design not frozen")
    rule=c["method"]
    if (rule["distance_bins_km"] != [0.4,1.2,2,5,10,20]
        or rule["allowed_original_csv_fields"] !=
           ["study","functional_group","taxon","source","site","experiment","fitness_type","lat_site","long_site"]):
        raise ValueError("Original metadata field or distance-bin contract drift")
    s=c["study"];st=s["Hargreaves_original_study_label"]; tax=s["original_organism"]
    origin=sorted(k for k in meta["origin_keys"] if k[0]==st and k[1]==tax)
    if len(origin)!=c["source"]["expected_original_source_population_ID_keys"]:
        raise ValueError("Original six strict-binomial source population labels changed")
    coords={}
    for key in origin:
        if key not in meta["gardens"]:
            raise ValueError("Source population not an exact study-species destination garden")
        pos=one_point(meta["gardens"][key])
        if pos is None:
            raise ValueError("Original garden/source coordinate is missing or ambiguous")
        coords[key]=pos
    if len(set(coords.values()))!=c["source"]["expected_unique_source_as_garden_proxy_points"]:
        raise ValueError("Original six distinct garden location points changed")
    # Source label may not encode published basin or quartz names. Never use
    # distances or latitude to invent assignments.
    codes={key:source_id_category(key[2],c) for key in origin}
    bcount=Counter(q["basin"] for q in codes.values())
    hcount=Counter(q["habitat"] for q in codes.values())
    declared=s["published_original_methods"]
    has_group_labels=(all(q["basin"] and q["habitat"] for q in codes.values())
         and all(bcount.get(b)==2 for b in declared["focal_drainage_basins"])
         and all(hcount.get(h)==3 for h in declared["quartz_habitat_groups"])
         and len({(q["basin"],q["habitat"]) for q in codes.values()})==6)
    pairdistance=[]
    within=[];between=[]
    for i,left in enumerate(origin):
        for right in origin[i+1:]:
            d=haversine(coords[left],coords[right])
            pairdistance.append(d)
            if has_group_labels:
                if codes[left]["basin"]==codes[right]["basin"]:
                    within.append(d)
                else:
                    between.append(d)
    if len(pairdistance)!=15:raise AssertionError("6 unique points must have 15 distinct unordered pairs")
    def summary(values):
        if not values:return None
        return {
          "count":len(values),
          "min_km":round(min(values),5),
          "median_km":round(percent_at(values,0.5),5),
          "max_km":round(max(values),5),
          "bins":dict(sorted(Counter(bucket(d,rule["distance_bins_km"]) for d in values).items()))
        }

    # Source-per-garden source design counts: experiment and fitness_type
    # qualify duplicate metadata cells, but do NOT increase source ID pairs.
    comparison_cells=[(key,sources) for key,sources in meta["source_groups"].items()
                      if key[0]==st and key[1]==tax and len(sources)>=2]
    if len(comparison_cells)!=c["source"]["expected_same_garden_multi_source_design_cells"]:
        raise ValueError("Previously frozen original same-garden comparison-cell count drift")
    gardenkeys={k[2] for k in origin}
    sourcekeys={k[2] for k in origin}
    directed=set()
    per_garden=defaultdict(set)
    unexpected_gardens=set()
    for cell,sources in comparison_cells:
        garden=cell[3]
        if garden not in gardenkeys:
            unexpected_gardens.add(garden)
        for source in sources:
            if source in sourcekeys and garden in gardenkeys:
                directed.add((source,garden))
                per_garden[garden].add(source)
    other_habitat_basin_edges=0
    if has_group_labels:
        info={k[2]:v for k,v in codes.items()}
        for source,garden in directed:
            if (info[source]["basin"]!=info[garden]["basin"]
                and info[source]["habitat"]!=info[garden]["habitat"]):
                other_habitat_basin_edges+=1
    result={
        "schema":"ttf_Ellis2007_SA_Knersvlakte_six_population_source_only_spatial_design_v01",
        "status":"SOURCE_ONLY_SIX_SITES_ORIGINAL_BARRIER_HABITAT_DESIGN_NOT_HISTORICAL_CLIMATE_VALIDATION",
        "original_paper_doi":s["original_doi"],
        "original_source_dataset_SHA256":c["source"]["original_Hargreaves_SHA256"],
        "study_original_binomial_taxa":1,
        "source_population_ID_keys":len(origin),
        "unique_garden_derived_proxy_location_points":len(set(coords.values())),
        "original_source_IDs_that_equal_same_study_garden_IDs":len(coords),
        "original_distinct_source_plus_garden_directed_ID_pairs":len(directed),
        "original_gardens_with_source_pair_metadata":len(per_garden),
        "per_garden_distinct_source_ID_count_distribution":
             dict(sorted(Counter(len(v) for v in per_garden.values()).items())),
        "original_source_metadata_comparison_cells_not_independent":len(comparison_cells),
        "original_unmatched_garden_metadata_label_count":len(unexpected_gardens),
        "original_study_basin_habitat_labels_directly_classifiable_from_ID_text":has_group_labels,
        "literal_basin_names_found_in_IDs":dict(sorted((str(k),v) for k,v in bcount.items())),
        "literal_quartz_habitat_names_found_in_IDs":dict(sorted((str(k),v) for k,v in hcount.items())),
        "all_15_unique_source_point_pair_distances_km":summary(pairdistance),
        "within_basin_3_pairs_if_literal_IDs_qualified":summary(within),
        "between_basin_12_pairs_if_literal_IDs_qualified":summary(between),
        "between_basin_different_quartz_habitat_edges_observed_if_codes_qualified":
             other_habitat_basin_edges if has_group_labels else None,
        "paper_provides_six_original_field_source_population_WGS84_points_as_exact_ID_table":False,
        "original_Hargreaves_destination_site_proxies_NOT_independent_source_GPS":True,
        "no_coordinate_or_source_ID_values_exported":True,
        "original_2000_adult_transplant_cohort_known_from_paper_methods":True,
        "source_ID_pair_geometry_is_not_15_independent_historical_exposures":True,
        "species_contrasts_are_edaphic_vs_basin_not_direct_LGM_selection_evidence":True,
        "original_fitness_or_genetic_outcomes_not_read_from_Hargreaves":True,
        "no_21ka_0BP_climate_predictor_or_fitness_effect_fitted":True,
        "Bontrager_terminal_44_of_50_gate_unchanged":True
    }
    return result

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--output",type=Path,required=True)
    args=p.parse_args()
    c=json.loads(CONTRACT.read_text(encoding="utf-8"))
    root=json.loads((ROOT/c["source"]["original_dataset_contract"]).read_text(encoding="utf-8"))
    entry=root["original"]["hargreaves"]
    raw=original_bytes(entry)
    if hashlib.sha256(raw).hexdigest()!=c["source"]["original_Hargreaves_SHA256"]:
        raise ValueError("Original source SHA256 drift")
    meta_contract=json.loads((ROOT/c["source"]["original_Hargreaves_source_only_parser_contract"]).read_text())
    meta=source_metadata(raw,"hargreaves",meta_contract)
    result=geometry(meta,c)
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(result,ensure_ascii=False,sort_keys=True,indent=2)+"\n")
    print("ELLIS2007_SIX_POPULATION_SOURCE_GEOGRAPHY="+json.dumps({
      "sourceID_count":result["source_population_ID_keys"],
      "source_garden_pairs":result["original_distinct_source_plus_garden_directed_ID_pairs"],
      "basin_habitat_labels_qualified":result["original_study_basin_habitat_labels_directly_classifiable_from_ID_text"],
      "pair_geometry":result["all_15_unique_source_point_pair_distances_km"],
      "within_basin":result["within_basin_3_pairs_if_literal_IDs_qualified"],
      "between_basin":result["between_basin_12_pairs_if_literal_IDs_qualified"],
      "no_LGM_fitness":result["no_21ka_0BP_climate_predictor_or_fitness_effect_fitted"]
    },sort_keys=True))
if __name__=="__main__":
    main()
