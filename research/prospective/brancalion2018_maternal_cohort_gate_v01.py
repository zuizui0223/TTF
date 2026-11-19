"""Brancalion 2018 parental-source cohort gate (source-only).

The three assay modalities share one harvest of the same maternal seed lots.
This is a source-unit/independence audit, not an analysis of fitness or SNPs.
"""
from __future__ import annotations
import argparse,json
from pathlib import Path

HERE=Path(__file__).resolve().parent
SHORTLIST=HERE/"hargreaves35_primary_article_recovery_shortlist_v01.json"
LABEL="Brancalionetal_EcolEvol_2018"
DOI="10.1002/ece3.4248"

PUBLISHED_METHODS={
    "forest_sites":["Caetetus Ecological Station","Carlos Botelho State Park",
                    "Ilha do Cardoso State Park"],
    "forest_populations":3,
    "mother_palms_per_population":10,
    "mother_palms_total":30,
    "harvesting_events":1,
    "seed_lots_shared_across_experiments":True,
    "assays":["common garden","seed reciprocal transplant",
              "seedling reciprocal transplant"],
}


def evaluate(shortlist):
    if shortlist.get("schema")!="ttf_Hargreaves35_original_primary_literature_recovery_shortlist_v01":
        raise ValueError("source-only shortlist changed")
    matching=[r for r in shortlist["record_list"]
              if r["original_Hargreaves_study_label"]==LABEL]
    if len(matching)!=1:raise ValueError("not exactly one Brancalion study label")
    r=matching[0]
    if (r["primary_article_doi"]!=DOI or
        r["source_only_candidate_same_garden_cells"]!=8 or
        r["distinct_original_source_ID_labels"]!=3 or
        r["distinct_original_garden_ID_labels"]!=3):
        raise ValueError("frozen Brancalion source geometry changed")
    if any(r[k] for k in ("primary_experiment_identity_certified",
                           "historical_climate_contrast_qualified",
                           "genetic_adaptation_plasticity_speciation_effect_opened")):
        raise ValueError("source-only study admission unexpectedly changed")
    methods=PUBLISHED_METHODS
    if (methods["forest_populations"]*methods["mother_palms_per_population"]
        !=methods["mother_palms_total"]):
        raise AssertionError("maternal source count inconsistent")
    if not methods["seed_lots_shared_across_experiments"]:
        raise AssertionError("source cohort unexpectedly separated")
    return {
       "schema":"ttf_Brancalion2018_single_maternal_cohort_source_gate_v01",
       "primary_paper_doi":DOI,
       "status":"SHARED_MATERNAL_COHORT_NO_INDEPENDENT_LGM_PANEL",
       "published_original_forest_sites":methods["forest_sites"],
       "published_original_source_populations":methods["forest_populations"],
       "published_mother_palms_per_population":methods["mother_palms_per_population"],
       "published_distinct_mother_palms_total":methods["mother_palms_total"],
       "published_seed_harvesting_events":methods["harvesting_events"],
       "experiment_modalities":methods["assays"],
       "experiment_modality_count":len(methods["assays"]),
       "modalities_share_original_maternal_seed_lots":True,
       "independent_maternal_source_cohorts_across_modalities":1,
       "hargreaves_same_garden_candidate_cells_not_independent":8,
       "hargreaves_original_source_ID_labels":3,
       "source_population_GPS_independently_verified":0,
       "experiment_cohort_IDs_independently_certified":False,
       "read_fitness_SNP_LGM_or_adaptation_response":False,
       "bontrager_44_of_50_stop_unchanged":True,
       "interpretation":"Same seed lots feed three distinct experimental modalities; no triple counting of independent source cohorts."
    }


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--output",type=Path,required=True)
    a=p.parse_args()
    r=evaluate(json.loads(SHORTLIST.read_text(encoding="utf-8")))
    a.output.parent.mkdir(parents=True,exist_ok=True)
    a.output.write_text(json.dumps(r,sort_keys=True,ensure_ascii=False,indent=2)+"\n",encoding="utf-8")
    print(json.dumps({k:r[k] for k in ("status","published_original_source_populations",
                                      "published_distinct_mother_palms_total",
                                      "independent_maternal_source_cohorts_across_modalities",
                                      "source_population_GPS_independently_verified")},sort_keys=True))


if __name__=="__main__":
    main()
