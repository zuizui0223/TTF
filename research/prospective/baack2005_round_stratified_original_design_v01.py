"""Baack & Stanton 2005 original three-round seedling transplant graph.

Source-only transcription of ORIGINAL Table 1 (Evolution 59:1936-1944).
Use site ID and cytotype only, not survival/growth, original GPS or LGM.
Three separate rounds must not be pooled as a complete 6x8 factorial.
"""
from __future__ import annotations
import argparse,json
from pathlib import Path

HERE=Path(__file__).resolve().parent
SHORTLIST=HERE/"hargreaves35_primary_article_recovery_shortlist_v01.json"
LABEL="BaackandStanton_Evolution_2005"
DOI="10.1111/j.0014-3820.2005.tb01063.x"
# Original published Table 1, in printed order. D=diploid, T=tetraploid.
ROUNDS={
  "Early":{"source":["D-LP","D-HPNG","T-HPWM","T-PM"],
           "garden":["D-LP","D-HPNG","T-HPWM","T-PM"]},
  "Middle":{"source":["D-LP","D-HPNG","T-HPWM","T-WP"],
            "garden":["D-HPNG","D-HPWW","T-P7","T-WP"]},
  "Late":{"source":["D-BP","D-HPNG","T-PM","T-WP"],
          "garden":["D-BP","D-HPNG","T-PM","T-WP"]},
}

def evaluate(shortlist):
    if shortlist.get("schema")!="ttf_Hargreaves35_original_primary_literature_recovery_shortlist_v01":
        raise ValueError("source-only original bibliography schema drift")
    candidates=[r for r in shortlist["record_list"]
                if r["original_Hargreaves_study_label"]==LABEL]
    if len(candidates)!=1:raise ValueError("original study label missing/duplicated")
    cand=candidates[0]
    if (cand["primary_article_doi"]!=DOI or
        cand["distinct_original_source_ID_labels"]!=6 or
        cand["distinct_original_garden_ID_labels"]!=8 or
        cand["source_only_candidate_same_garden_cells"]!=8):
        raise ValueError("original metadata/Hargreaves comparison scope changed")
    if any(cand[k] for k in ("primary_experiment_identity_certified",
                              "historical_climate_contrast_qualified",
                              "genetic_adaptation_plasticity_speciation_effect_opened")):
        raise ValueError("outcome-blind qualification unexpectedly changed")
    all_rows=[]
    per_round=[]
    for round_name,design in ROUNDS.items():
        sources=design["source"];gardens=design["garden"]
        if len(sources)!=4 or len(gardens)!=4 or len(set(sources))!=4 or len(set(gardens))!=4:
            raise ValueError("published round design changed")
        if sum(x.startswith("D-") for x in sources)!=2 or sum(x.startswith("D-") for x in gardens)!=2:
            raise ValueError("cytotype balance changed")
        edges=[(round_name,s,g) for s in sources for g in gardens]
        all_rows+=edges
        per_round.append({
           "round":round_name,
           "source_populations":len(sources),
           "garden_populations":len(gardens),
           "directed_transplant_arms":len(edges),
           "home_arms":sum(s==g for _,s,g in edges),
           "cross_cytotype_environment_arms":sum(s[0]!=g[0] for _,s,g in edges),
           "sources_omitted_from_some_other_rounds":True,
        })
    sources={s for _,s,_ in all_rows}
    gardens={g for _,_,g in all_rows}
    unique_pairs={(s,g) for _,s,g in all_rows}
    if (len(all_rows),len(unique_pairs),len(sources),len(gardens),
        sum(x["home_arms"] for x in per_round))!=(48,38,6,8,10):
        raise AssertionError("original Table 1 reconstruction drift")
    return {
      "schema":"ttf_Baack2005_original_seedling_round_stratified_design_v01",
      "original_article_doi":DOI,
      "status":"ROUND_SPECIFIC_EXPERIMENT_DESIGN_NOT_COMPLETE_6x8_MATRIX_OR_LGM_PANEL",
      "source_only_hargreaves_candidate_cells":8,
      "original_study_unique_populations_reported":10,
      "original_seedling_source_populations":len(sources),
      "original_seedling_garden_populations":len(gardens),
      "ploidies_compared":["D","T"],
      "seedling_experiment_rounds":per_round,
      "round_specific_directed_transplant_arms":len(all_rows),
      "distinct_source_garden_pairs_after_ignoring_round":len(unique_pairs),
      "home_arms_with_round_preserved":10,
      "cross_cytotype_environment_arms_with_round_preserved":sum(s[0]!=g[0] for _,s,g in all_rows),
      "complete_6x8_factorial_present":False,
      "source_GPS_originally_certified":0,
      "Hargreaves_cells_not_directly_equivalent_to_round_stratified_arms":True,
      "original_seed_experiment_and_controlled_chamber_not_counted_as_independent_studies":True,
      "ploidy_and_local_environment_must_not_be_misattributed_to_LGM":True,
      "fitness_genetic_adaptation_plasticity_speciation_or_LGM_responses_opened":False,
      "old_bontrager_44_of_50_stop_unchanged":True,
      "original_table1_round_populations":ROUNDS,
    }

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--output",required=True,type=Path)
    args=p.parse_args()
    d=json.loads(SHORTLIST.read_text("utf-8"))
    r=evaluate(d)
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(r,sort_keys=True,indent=2)+"\n",encoding="utf-8")
    print(json.dumps({k:v for k,v in r.items()
                      if k not in ("original_table1_round_populations","seedling_experiment_rounds")},
                     sort_keys=True))
if __name__=="__main__":
    main()
