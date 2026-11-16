"""Tomiolo source-only garden-derived points × published appendix site geography.

Uses ONLY Hargreaves original exact (study,taxon,source) strings plus same
exact destination garden ID WGS84, no response or fitness fields. Report
anonymized 2×3 distances to published whole-arcminute site points M, SA, A.
Distances are exploratory site concordance, NOT source population coordinates
or independent historical environmental exposure measurements.
"""
from __future__ import annotations
import argparse,hashlib,json,math
from collections import defaultdict
from pathlib import Path
from bontrager_hargreaves_finegrain_overlap_v01 import original_bytes,haversine
from hargreaves_source_as_garden_georecovery_v01 import source_metadata,one_point

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
CONTRACT=HERE/"tomiolo2015_geo_proxy_site_matrix_contract_v01.json"

def proxy_geography(metadata:dict,c:dict)->dict:
    s=c["source"];geometry=c["from_previously_SHA_frozen_original_metadata"]
    r=c["matching_rules"]
    if c["status"]!="FROZEN_BEFORE_ORIGINAL_HARGREAVES_SITE_DISTANCE_MATRIX_EXECUTION":
        raise ValueError("No frozen Tomiolo original geography design")
    if (r["descriptive_distance_bands_km"]!=[0.5,1,2,5,10,50]
        or r["source_population_ID_count_must_equal_6"] is not True
        or not r["report_full_2_by_3_Haversine_distance_matrix_km_not_nearest_site_ID_assignment"]):
        raise ValueError("Original source-only no-retune rules drifted")
    study=s["Hargreaves_original_study_label"]
    keys={k for k in metadata["origin_keys"] if k[0]==study}
    if len(keys)!=geometry["original_study_plant_sourceID_keys_across_species"]:
        raise ValueError("Tomiolo original source population ID key count changed")
    pts={}
    by_taxon=defaultdict(set)
    for key in keys:
        if key not in metadata["gardens"]:
            raise ValueError("Original source population ID not also a same-text garden")
        candidate=one_point(metadata["gardens"][key])
        if candidate is None:
            raise ValueError("Source-as-garden geography is missing or ambiguous")
        pts[key]=candidate
        by_taxon[key[1]].add(candidate)
    unique=sorted(set(pts.values()),key=lambda pt:(pt[1],pt[0]))
    if (len(unique)!=geometry["original_study_unique_garden_proxy_coordinate_points"]
        or len(pts)!=geometry["original_study_distinct_sourceIDs_with_garden_proxy"]):
        raise ValueError("Previously frozen source-only garden-coordinate geometry drifted")
    names=["M","SA","A"]
    pub=c["published_site_locations"]
    if set(pub)!=set(names) or len(names)!=3:
        raise ValueError("Published original climate site table changed")
    distances={}
    for i,proxy in enumerate(unique,1):
        item={}
        for name in names:
            loc=pub[name]
            xy=(loc["lon"],loc["lat"])
            if not all(math.isfinite(float(v)) for v in xy):
                raise ValueError("Original published site coordinate invalid")
            item[name]=round(haversine(proxy,xy),4)
        distances[f"anonymous_proxy_point_{i}"]=item
    # Produce full distance matrix; never pick a matching M/SA label for a
    # source ID or change any source ID based on geographic closeness.
    return {
       "schema":"ttf_Tomiolo2015_two_garden_derived_sites_three_published_appendix_sites_geo_matrix_v01",
       "status":"EXPLORATORY_TWO_GEOGRAPHIC_POINTS_VS_PUBLISHED_SITES_NOT_GPS_SOURCE_PROOF",
       "original_paper_doi":s["paper_doi"],
       "original_Hargreaves_study_label":study,
       "original_study_source_population_ID_keys":len(keys),
       "original_study_unique_garden_derived_site_points":len(unique),
       "original_study_binomial_species_labels":len(by_taxon),
       "strict_binomial_species_with_two_distinct_garden_proxy_locations":
           sum(len(v)>=2 for v in by_taxon.values()),
       "published_author_site_labels":["M","SA","A"],
       "original_paper_field_origin_sites":["M","SA"],
       "A_is_reference_aridity_gradient_site_not_automatically_an_experiment_source":True,
       "anonymous_proxy_point_to_published_site_distance_km":distances,
       "nominal_published_coordinate_precision":"whole arcminute",
       "study_site_code_to_source_population_ID_directly_certified":False,
       "Hargreaves_garden_proxy_coords_exported":False,
       "source_population_GPS_directly_verified":False,
       "independent_primary_experiment_cohort_verified":False,
       "no_fitness_phenotype_genotype_LGM_model_opened":True,
       "old_Bontrager_44of50_source_gate_unchanged":True
    }

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--output",type=Path,required=True)
    args=p.parse_args()
    c=json.loads(CONTRACT.read_text(encoding="utf-8"))
    base=json.loads((ROOT/c["source"]["Hargreaves_original_source_metadata_contract"]).read_text(encoding="utf-8"))
    entry=base["original"]["hargreaves"]
    data=original_bytes(entry)
    if hashlib.sha256(data).hexdigest()!=c["source"]["Hargreaves_exact_sha256"]:
        raise ValueError("Authentic original source archive changed")
    source_rule=json.loads((ROOT/c["source"]["source_as_garden_metadata_contract"]).read_text())
    meta=source_metadata(data,"hargreaves",source_rule)
    result=proxy_geography(meta,c)
    result["source_original_SHA256_verified"]=c["source"]["Hargreaves_exact_sha256"]
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(result,indent=2,ensure_ascii=False,sort_keys=True)+"\n",encoding="utf-8")
    print("TOMIOLO2015_TWO_SITE_SOURCE_GEO="+json.dumps({
      "original_source_IDs":result["original_study_source_population_ID_keys"],
      "unique_garden_proxy_points":result["original_study_unique_garden_derived_site_points"],
      "species_with_two_proxy_sites":result["strict_binomial_species_with_two_distinct_garden_proxy_locations"],
      "anonymous_distance_matrix_km":result["anonymous_proxy_point_to_published_site_distance_km"],
      "population_GPS_certified":result["source_population_GPS_directly_verified"]
    },sort_keys=True))

if __name__=="__main__":
    main()
