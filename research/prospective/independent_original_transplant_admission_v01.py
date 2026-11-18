"""Response-blind eligibility audit for independent original transplant panels.

Uses only a frozen bibliographic shortlist, never outcome values. A DOI,
source label, or source/garden cell count does not certify independent
original experiments or original source GPS.
"""
from __future__ import annotations
import argparse,hashlib,json
from pathlib import Path

HERE=Path(__file__).resolve().parent
SOURCE=HERE/"hargreaves35_primary_article_recovery_shortlist_v01.json"
SOURCE_BLOB_SHA="d21413b4900ee7fef84439495b26fc199864c315"
PRIORITY={
 "Brancalionetal_EcolEvol_2018":"PRIMARY_POPULATION_GPS_AND_COHORT_ID",
 "BaackandStanton_Evolution_2005":"PRIMARY_POPULATION_GPS_PLOIDY_AND_COHORT_ID",
 "Kalskeetal_JEcol_2012":"ORIGINAL_DMS_COORDINATE_VALIDATION",
 "Ellisetal_NewPhytol_2007":"ORIGINAL_MISSING_SIX_ARMS_AND_SOURCE_GPS",
 "Raabovaetal_BiolConservation_2007":"ORIGINAL_ADULT_COHORT_ID_AND_GPS",
 "Tomioloetal_Ecology_2015":"ONLY_TWO_REPRESENTATIVE_SITES_NOT_INDEPENDENT_GPS",
 "Meloetal_NewPhytol_2014":"PRIMARY_SOURCE_GPS_AND_EXPERIMENT_ID",
 "Volisetal_Oecologia_2002":"PRIMARY_SOURCE_GPS_AND_EXPERIMENT_ID",
 "Aranyetal_PlantEcol_2009":"PRIMARY_SOURCE_GPS_AND_EXPERIMENT_ID",
 "Kindelletal_JEcol_1996":"PRIMARY_SOURCE_GPS_AND_EXPERIMENT_ID"
}

def evaluate(source):
    if source.get("schema")!="ttf_Hargreaves35_original_primary_literature_recovery_shortlist_v01":
        raise ValueError("source shortlist schema changed")
    rows=source["record_list"]
    if len(rows)!=10 or {r["original_Hargreaves_study_label"] for r in rows}!=set(PRIORITY):
        raise ValueError("source study set drift")
    if source["summary"]["study_labels_with_original_source_lat_lon_verified"]!=0 or source["summary"]["independent_experiments_certified"]!=0:
        raise ValueError("unexpected source qualification; re-audit required")
    for r in rows:
        if not r["first_author_year_journal_title_bibliographic_match"]:
            raise ValueError("unverified DOI/title")
        if r["primary_experiment_identity_certified"] or r["historical_climate_contrast_qualified"] or r["genetic_adaptation_plasticity_speciation_effect_opened"]:
            raise ValueError("frozen source-only boundary changed")
        if r["source_only_candidate_same_garden_cells"]<0:
            raise ValueError("negative comparison count")
    return {
      "schema":"ttf_IndependentTransplant_original_source_admission_matrix_v01",
      "status":"NO_NEW_INDEPENDENT_HISTORICAL_FITNESS_PANEL_QUALIFIED",
      "source_shortlist_git_blob_sha":SOURCE_BLOB_SHA,
      "bibliographic_labels_checked":len(rows),
      "candidate_cells_not_independent_replicates":sum(r["source_only_candidate_same_garden_cells"] for r in rows),
      "certified_original_source_population_GPS_studies":0,
      "certified_independent_original_experiments":0,
      "historical_source_contrast_qualified_studies":0,
      "first_next_action":"ORIGINAL_SOURCE_GPS_AND_EXPERIMENT_ID_RECOVERY_NOT_FITNESS_OPENING",
      "priority_requests":[
        {"study_label":r["original_Hargreaves_study_label"],
         "doi":r["primary_article_doi"],
         "metadata_needed":PRIORITY[r["original_Hargreaves_study_label"]],
         "status":"SOURCE_ONLY_HOLD"}
        for r in rows
      ],
      "prior_bontrager_44_of_50_stop_unchanged":True,
      "fitness_genetics_plasticity_speciation_or_LGM_response_opened":False
    }

def main():
    p=argparse.ArgumentParser();p.add_argument("--output",type=Path,required=True)
    a=p.parse_args()
    r=evaluate(json.loads(SOURCE.read_text(encoding="utf-8")))
    a.output.parent.mkdir(parents=True,exist_ok=True)
    a.output.write_text(json.dumps(r,indent=2,sort_keys=True)+"\n",encoding="utf-8")
    print(json.dumps({k:v for k,v in r.items() if k!="priority_requests"},sort_keys=True))

if __name__=="__main__":main()
