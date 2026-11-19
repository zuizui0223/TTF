"""Publication-versus-experiment-unit deduplication from original-source evidence only.

One author-confirmed adult transplant started in October 2004 has a 2007 and
2014 report. Raabova's 2003 and 2004 reciprocal *seed sowing* experiments are
separate original interventions and must never be merged into that adult
transplant simply because they appear in the same 2007 paper.

No original fitness, genetic, species-divergence or historical-climate
measurements are read or analyzed. This is a nonresponse lineage example.
"""
from __future__ import annotations
import argparse
import json
from collections import Counter
from pathlib import Path

HERE = Path(__file__).resolve().parent
EVIDENCE = HERE / "raabova2007_pankova2014_author_confirmed_adult_cohort_v02.json"

def deduplicate_cohorts(evidence: dict) -> dict:
    if evidence.get("status") != "AUTHOR_EXPLICIT_SHARED_ADULT_2004_TRANSPLANT_EXPERIMENT_INDIVIDUAL_ID_UNVERIFIED":
        raise ValueError("No publisher-explicit same-adult-experiment source evidence")
    src = evidence["evidence"]
    if (src["reference_number"] != 21
        or src["reference_doi"] != "10.1016/j.biocon.2007.07.007"
        or src["followup_doi"] != "10.1371/journal.pone.0093967"
        or not src["authors_explicitly_identify_same_experiment"]
        or src["short_source_quote"] != "same reciprocal transplant experiments after two years"):
        raise ValueError("Author bibliography link changed; do not assume same experiment")
    c = evidence["cohort"]
    if (c["field_transplant"] != "October 2004"
        or c["original_six_origin_and_garden_codes"] != ["K1","K2","K3","S1","S2","S3"]
        or c["plants_per_origin_garden"] != 25
        or not c["separate_year2003_year2004_seed_sowing_experiments_not_merged"]):
        raise ValueError("Original adult-versus-seed intervention boundaries changed")
    if (evidence["deduplication"]["shared_adult_cohort_count"] != 1
        or evidence["deduplication"]["incremental_independent_adult_cohorts_from_followup"] != 0
        or evidence["deduplication"]["raw_individual_ID_crosspaper_join_verified"]):
        raise ValueError("Prior experimental cohort qualification drifted")
    paper2007, paper2014 = src["reference_doi"],src["followup_doi"]
    records = [
        {"paper_doi":paper2007, "reported_component":"reciprocal_seed_sowing_2003",
         "underlying_experiment_id":"Raabova2007_seed_sowing_2003", "time_reported":"2003–2006"},
        {"paper_doi":paper2007, "reported_component":"reciprocal_seed_sowing_2004",
         "underlying_experiment_id":"Raabova2007_seed_sowing_2004", "time_reported":"2004–2006"},
        {"paper_doi":paper2007, "reported_component":"adult_field_transplant_2004",
         "underlying_experiment_id":c["key"],"time_reported":"2005–2006"},
        {"paper_doi":paper2014,"reported_component":"adult_field_transplant_followup_2009",
         "underlying_experiment_id":c["key"],"time_reported":"2009"},
    ]
    report_components = {(x["paper_doi"],x["reported_component"]) for x in records}
    unique = {x["underlying_experiment_id"] for x in records}
    adult = [x for x in records if x["underlying_experiment_id"]==c["key"]]
    if len(report_components)!=4 or len(unique)!=3 or len(adult)!=2:
        raise AssertionError("Experiment/report lineage graph internal inconsistency")
    if len({x["paper_doi"] for x in adult})!=2:
        raise AssertionError("Original and follow-up papers not both present")
    if len({x["underlying_experiment_id"] for x in adult})!=1:
        raise AssertionError("Same adult field cohort falsely split")
    return {
       "schema":"ttf_author_confirmed_2004_experiment_component_dedup_source_only_v01",
       "status":"ONE_PHYSICAL_ADULT_EXPERIMENT_TWO_PAPERS_SEED_SOWING_COMPONENTS_DISTINCT",
       "report_components":records,
       "original_unique_article_DOIs":len(set(x["paper_doi"] for x in records)),
       "article_by_experimental_component_records":len(records),
       "distinct_underlying_experiment_cohorts":len(unique),
       "adult_transplant_report_papers":len(adult),
       "adult_transplant_underlying_cohorts":len({x["underlying_experiment_id"] for x in adult}),
       "independent_adult_experiments_added_by_later_2014_paper":0,
       "original_2007_seed_sowing_intervention_cohorts":2,
       "same_species_or_region_does_not_imply_same_experiment":True,
       "author_statement_direct_not_only_author_year_or_DOI_similarity":True,
       "original_individual_tag_IDs_directly_matched":False,
       "raw_fitness_or_genetic_observations_read":False,
       "previous_44of50_LGM_adaptation_source_gate_unchanged":True,
       "does_not_certify_three_statistically_independent_species_or_populations":True,
    }

def main() -> None:
    p=argparse.ArgumentParser()
    p.add_argument("--output", type=Path,required=True)
    args=p.parse_args()
    evidence=json.loads(EVIDENCE.read_text(encoding="utf-8"))
    output=deduplicate_cohorts(evidence)
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(output,ensure_ascii=False,indent=2,sort_keys=True)+"\n",encoding="utf-8")
    print("TTF_EXPERIMENT_REPORT_DEDUP=",json.dumps({
        "article_DOIs":output["original_unique_article_DOIs"],
        "report_components":output["article_by_experimental_component_records"],
        "physical_experiment_cohorts":output["distinct_underlying_experiment_cohorts"],
        "adult_followup_increment":output["independent_adult_experiments_added_by_later_2014_paper"],
        "status":output["status"],
    },sort_keys=True))

if __name__=="__main__":
    main()
