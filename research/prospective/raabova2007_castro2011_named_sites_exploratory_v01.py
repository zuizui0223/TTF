"""Exploratory correctly-referenced Raabova2007 S1-S3 vs Castro2011 site QA.

A *new, descriptive source-only* check. 2010 Sudova table does NOT cite 2007.
Match 2007 original S1, S2, S3 labels by literal code inside original study,
not by nearest coordinates, and compare to 2011 diploid site NAMES.
Never use to rescue a failed confirmatory geographic qualification.
"""
from __future__ import annotations
import argparse,hashlib,json
from pathlib import Path

from bontrager_hargreaves_finegrain_overlap_v01 import original_bytes,haversine
from hargreaves_source_as_garden_georecovery_v01 import source_metadata,one_point
from raabova_2007_2008_original_ID_metadata_v01 import classify_code
from raabova_sudova_crosspaper_sites_v01 import dms

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
CONTRACT=HERE/"raabova2007_castro2011_named_sites_exploratory_contract_v01.json"

def match_by_verified_primary_population_codes(obs:dict,c:dict,erratum:dict)->dict:
    if c["status"]!="EXPLORATORY_SOURCE_ONLY_FROZEN_BEFORE_NEW_PUBLISHED_SITE_COMPARISON":
        raise ValueError("New exploratory source-only contract not frozen")
    primary=c["primary_original"]
    other=c["corroborating_external_publication"]
    match=c["matching_rule"]
    err=erratum["author_publication_referent"]
    if (err["referenced_primary_DOI"]!="10.1007/s00442-008-1156-1"
        or err["erroneously_assigned_2007_DOI"]!=primary["paper_doi"]
        or other["paper_doi"]!="10.1007/s10682-010-9439-5"):
        raise ValueError("Wrong literature identity or invalid original referent")
    if (match["descriptive_cutoff_km"]!=1
        or not match["no_fuzzy_alias_or_nearest_site_number_matching"]
        or not match["descriptive_1km_cutoff_only_not_a_confirmatory_acceptance_rule"]):
        raise ValueError("Exploratory distance/matching contract changed")
    study=primary["Hargreaves_study_label"]
    taxon=primary["taxon"]
    allowed=primary["published_2007_population_codes"]
    origin={k for k in obs["origin_keys"] if k[0]==study and k[1]==taxon}
    if len(origin)!=6:
        raise ValueError("Original Raabova2007 six source population identifiers drifted")
    code_to_coord={}
    for key in origin:
        code=classify_code(key[2],allowed)
        if code not in allowed or code in code_to_coord:
            raise ValueError("Original Raabova2007 population codes not unique and literal")
        if key not in obs["gardens"]:
            raise ValueError("No identical same-study garden identifier for original source code")
        point=one_point(obs["gardens"][key])
        if point is None:
            raise ValueError("Ambiguous/missing origin-as-garden coordinate")
        code_to_coord[code]=point
    if set(code_to_coord)!=set(allowed):
        raise ValueError("Original six-population K/S codes incomplete")
    published=other["Table1_diploid_populations"]
    if [x["primary_2007_population_code"] for x in published]!=["S1","S2","S3"]:
        raise ValueError("Externally published named diploid localities changed or were reordered")
    if [x["site_name"] for x in published]!=[primary["S1_name"],primary["S2_name"],primary["S3_name"]]:
        raise ValueError("Primary 2007 and 2011 site-name evidence drifted")
    output=[]
    for record in published:
        code=record["primary_2007_population_code"]
        external=(dms(record["longitude_DMS"],"lon"),
                  dms(record["latitude_DMS"],"lat"))
        dist=haversine(code_to_coord[code],external)
        output.append({
            "verified_Raabova2007_source_population_code":code,
            "named_locality_shared_in_2007_and_2011":record["site_name"],
            "source_proxy_to_2011_pub_coordinate_km":round(dist,6),
            "descriptive_distance_le_original_1km_cutoff":dist<=match["descriptive_cutoff_km"],
            "original_Hargreaves_source_ID_and_coordinate_not_exported":True,
        })
    return {
        "schema":"ttf_Raabova2007_Castro2011_diploid_sites_exploratory_name_matched_geodistance_v01",
        "status":"EXPLORATORY_VALID_PAPER_AND_NAME_CROSSWALK_NOT_INDEPENDENT_GPS_CERTIFICATION",
        "target_Raabova_2007_DOI":primary["paper_doi"],
        "corroborating_Castro_2011_DOI":other["paper_doi"],
        "corroborating_Sudova2010_2008_wrong_reference_EXCLUDED":True,
        "Hargreaves_2007_all_six_literal_population_codes_present":True,
        "published_named_diploid_localities_compared":len(output),
        "original_1km_distance_threshold_unchanged":True,
        "exploratory_published_name_matched_sites":output,
        "within_1km_descriptive_site_count":sum(x["descriptive_distance_le_original_1km_cutoff"] for x in output),
        "scientific_scope":"One exploratory source-only matched-name site check; these published GPS may share source georeferencing or represent different microplots.",
        "independently_verified_raw_field_GPS_coordinates":0,
        "independent_primary_experiments_certified":0,
        "no_fitness_genetic_LGM_adaptation_plasticity_speciation_outcomes_fitted":True,
        "preexisting_frozen_44_of_50_adaptation_stop_unchanged":True,
        "old_invalid_1_of_4_external_test_remains_invalid":True,
    }

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--output",type=Path,required=True)
    a=p.parse_args()
    c=json.loads(CONTRACT.read_text(encoding="utf-8"))
    errpath=ROOT/c["original_source"]["correct_referent_erratum"]
    err=json.loads(errpath.read_text(encoding="utf-8"))
    if err["status"]!="TERMINAL_INVALID_EXTERNAL_VALIDATION_WRONG_PRIMARY_PAPER_REFERENT":
        raise ValueError("Prior incorrect scientific reference not closed")
    original_contract=json.loads((ROOT/c["original_source"]["original_source_contract"]).read_text(encoding="utf-8"))
    entry=original_contract["original"]["hargreaves"]
    raw=original_bytes(entry)
    if hashlib.sha256(raw).hexdigest()!=c["original_source"]["original_Hargreaves_SHA256"]:
        raise ValueError("Original Hargreaves source input changed")
    asgard_rule=json.loads((HERE/"hargreaves_source_as_garden_georecovery_contract_v01.json").read_text())
    obs=source_metadata(raw,"hargreaves",asgard_rule)
    r=match_by_verified_primary_population_codes(obs,c,err)
    r["original_source_SHA256"]=c["original_source"]["original_Hargreaves_SHA256"]
    a.output.parent.mkdir(parents=True,exist_ok=True)
    a.output.write_text(json.dumps(r,indent=2,ensure_ascii=False,sort_keys=True)+"\n",encoding="utf-8")
    print("RAABOVA2007_CASTRO2011_CORRECT_SITE_EXPLORATORY=",json.dumps({
        "status":r["status"],
        "sites":r["exploratory_published_name_matched_sites"],
        "descriptive_within_1km":r["within_1km_descriptive_site_count"]
    },ensure_ascii=False,sort_keys=True))

if __name__=="__main__":main()
