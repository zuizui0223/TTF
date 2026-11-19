"""Two remaining 10-shortlist original study designs: Volis 2002 and Kindell 1996.

SOURCE-ONLY. No source GPS, LGM, fitness, genetics, or original individual
outcome data. Volis Methods: author-hosted 2002 Oecologia original PDF.
Kindell: published 1996 abstract only; full original arm schedule unverified.
"""
from __future__ import annotations
import argparse,json
from pathlib import Path

HERE=Path(__file__).resolve().parent
SHORTLIST=HERE/"hargreaves35_primary_article_recovery_shortlist_v01.json"
EXPECTED={
    "Volisetal_Oecologia_2002":{
        "doi":"10.1007/s00442-002-0999-0",
        "source_labels_in_hargreaves":2,
        "garden_labels_in_hargreaves":2,
        "candidate_cells":10,
        "original_source_populations":2,
        "original_field_gardens":2,
        "source_data_level":"AUTHOR_HOSTED_ORIGINAL_OECOLOGIA_PDF_METHODS",
        "original_methods_url":"https://tzin.bgu.ac.il/life/Faculty/Volis/pdf/pdf4.pdf",
        "distinct_site_contexts":["Beit Guvrin Mediterranean population",
                                  "Negev wadi desert population"],
        "cleared_vegetation_maternal_parents_per_origin":3,
        "intact_vegetation_maternal_parents_per_origin":10,
        "distinct_maternal_subexperiment_pools":2,
        "separate_original_allozyme_10_population_paper_is_not_this_two_site_transplant":True,
        "original_age_treatments":None,
        "vegetation_treatments":"intact versus cleared; planting density and mixture varied",
    },
    "Kindelletal_JEcol_1996":{
        "doi":"10.2307/2261336",
        "source_labels_in_hargreaves":2,
        "garden_labels_in_hargreaves":6,
        "candidate_cells":6,
        "original_source_populations":2,
        "original_field_gardens":2,
        "source_data_level":"PUBLISHED_ABSTRACT_NOT_VERIFIED_FULL_ORIGINAL_METHODS",
        "original_methods_url":"https://www.jstor.org/stable/2261336",
        "distinct_site_contexts":["original perennial grass population A (unlocalized)",
                                  "original perennial grass population B (unlocalized)"],
        "cleared_vegetation_maternal_parents_per_origin":None,
        "intact_vegetation_maternal_parents_per_origin":None,
        "distinct_maternal_subexperiment_pools":None,
        "separate_original_allozyme_10_population_paper_is_not_this_two_site_transplant":None,
        "original_age_treatments":3,
        "vegetation_treatments":"undisturbed versus removed vegetation",
    },
}

def evaluate(d):
    if d.get("schema")!="ttf_Hargreaves35_original_primary_literature_recovery_shortlist_v01":
        raise ValueError("Frozen bibliography schema changed")
    rows=d["record_list"]
    if len(rows)!=10 or len({x["original_Hargreaves_study_label"] for x in rows})!=10:
        raise ValueError("Frozen 10-study admission panel changed")
    items=[]
    for label,gate in EXPECTED.items():
        matches=[x for x in rows if x["original_Hargreaves_study_label"]==label]
        if len(matches)!=1:raise ValueError("Study label missing/duplicated")
        r=matches[0]
        if (r["primary_article_doi"]!=gate["doi"]
            or r["distinct_original_source_ID_labels"]!=gate["source_labels_in_hargreaves"]
            or r["distinct_original_garden_ID_labels"]!=gate["garden_labels_in_hargreaves"]
            or r["source_only_candidate_same_garden_cells"]!=gate["candidate_cells"]):
            raise ValueError("Original-vs-Hargreaves study grain mismatch")
        if any(r[k] for k in ("primary_experiment_identity_certified",
                              "historical_climate_contrast_qualified",
                              "genetic_adaptation_plasticity_speciation_effect_opened")):
            raise ValueError("Frozen biological qualification changed")
        assert gate["original_source_populations"]==2
        if label.startswith("Volis"):
            if gate["cleared_vegetation_maternal_parents_per_origin"]!=3:
                raise AssertionError("Cleared vegetation original seed pool drift")
            if gate["intact_vegetation_maternal_parents_per_origin"]!=10:
                raise AssertionError("Intact vegetation original seed pool drift")
        else:
            if gate["original_age_treatments"]!=3 or gate["garden_labels_in_hargreaves"]!=6:
                raise AssertionError("Age or garden labels drift")
        items.append({
            "study_label":label,"doi":gate["doi"],
            "evidence_resolution":gate["source_data_level"],
            "original_evidence_url":gate["original_methods_url"],
            "source_habitat_contexts":gate["distinct_site_contexts"],
            "original_source_populations":gate["original_source_populations"],
            "original_field_garden_sites":gate["original_field_gardens"],
            "Hargreaves_source_labels":gate["source_labels_in_hargreaves"],
            "Hargreaves_garden_labels_NOT_certified_as_distinct_field_geographies":
                gate["garden_labels_in_hargreaves"],
            "Hargreaves_candidate_cells_NOT_independent_original_studies":gate["candidate_cells"],
            "published_source_only_design_traits":{
                "original_age_treatments":gate["original_age_treatments"],
                "vegetation_treatments":gate["vegetation_treatments"],
                "cleared_vegetation_seed_pool_parents_per_origin":
                    gate["cleared_vegetation_maternal_parents_per_origin"],
                "intact_vegetation_seedling_parents_per_origin":
                    gate["intact_vegetation_maternal_parents_per_origin"],
            },
            "complete_primary_arm_schedule_authenticated":False,
            "original_population_source_GPS_certified":0,
            "independent_historical_climate_contrast_qualified":False,
            "fitness_genetic_adaptation_plasticity_speciation_LGM_response_opened":False,
            "status":"SOURCE_ONLY_HOLD_NOT_NEW_INDEPENDENT_LGM_PANEL",
        })
    return {
        "schema":"ttf_Volis2002_Kindell1996_original_method_scope_closure_v01",
        "status":"TEN_SHORTLIST_PAPERS_SOURCE_ONLY_REVIEWED_ZERO_HISTORICAL_FITNESS_PANELS_QUALIFIED",
        "original_shortlist_label_count":10,
        "newly_scoped_last_two_original_studies":2,
        "studies":items,
        "full_original_methods_verified_for_kindell":False,
        "two_original_source_populations_each":True,
        "existing_shortlist_verified_original_source_GPS_studies":0,
        "existing_shortlist_independent_original_experiments_certified":0,
        "historical_source_contrast_qualified":0,
        "independent_10_paper_panel_qualified":False,
        "do_not_count_treatments_or_separate_related_papers_as_extra_origins":True,
        "no_fitness_genetics_LGM_adaptation_plasticity_speciation_response_opened":True,
        "prior_Bontrager_44of50_stop_unchanged":True,
    }

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--output",required=True,type=Path)
    a=p.parse_args()
    r=evaluate(json.loads(SHORTLIST.read_text(encoding="utf-8")))
    a.output.parent.mkdir(parents=True,exist_ok=True)
    a.output.write_text(json.dumps(r,indent=2,sort_keys=True)+"\n",encoding="utf-8")
    print(json.dumps({k:v for k,v in r.items() if k!="studies"},sort_keys=True))

if __name__=="__main__":
    main()
