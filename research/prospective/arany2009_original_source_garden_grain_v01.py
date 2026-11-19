"""Arany et al. 2009 original source/garden provenance and replicate-grain audit.

Original Plant Ecology DOI 10.1007/s11258-008-9530-y; Methods only.
Published source population labels = 3 dune + 2 inland, but neither
a mother plant nor a descendant rosette is an extra source-GPS unit.
Inland enclosure is a *nearby common garden*, not one of the original
inland seed collection sites. No biological response rows are read.
"""
from __future__ import annotations
import argparse
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
SHORTLIST = HERE / "hargreaves35_primary_article_recovery_shortlist_v01.json"
LABEL = "Aranyetal_PlantEcol_2009"
DOI = "10.1007/s11258-008-9530-y"
METHODS = {
    "dune_original_seed_populations": 3,
    "inland_original_seed_populations": 2,
    "mother_plants_sampled_per_population": 10,
    "enclosures": 2,
    "descendant_rosettes_per_mother_per_enclosure": 8,
    "distinct_original_source_habitat_classes": 2,
    "inland_garden_is_original_seed_collection_site": False,
    "dune_enclosure_near_original_dune_population_3": True,
    "inland_enclosure_near_inland_population_1_but_not_at_original_inland_sites": True,
    "dune_population_1_to_2_reported_distance_m": 5,
    "dune_pop3_to_other_dune_approx_distance_m": 500,
    "inland_1_to_2_approx_distance_km": 8,
    "dune_to_inland_min_approx_distance_km": 6,
}


def evaluate(data):
    if data.get("schema") != "ttf_Hargreaves35_original_primary_literature_recovery_shortlist_v01":
        raise ValueError("Original source registry schema changed")
    found = [x for x in data["record_list"]
             if x["original_Hargreaves_study_label"] == LABEL]
    if len(found) != 1:
        raise ValueError("Arany original label missing or duplicated")
    r = found[0]
    if (r["primary_article_doi"] != DOI
        or r["distinct_original_source_ID_labels"] != 5
        or r["distinct_original_garden_ID_labels"] != 2
        or r["source_only_candidate_same_garden_cells"] != 6):
        raise ValueError("Hargreaves source/garden label drift")
    if any(r[k] for k in (
        "primary_experiment_identity_certified",
        "historical_climate_contrast_qualified",
        "genetic_adaptation_plasticity_speciation_effect_opened",
    )):
        raise ValueError("Original source-only qualification changed")
    m = METHODS
    sources = m["dune_original_seed_populations"] + m["inland_original_seed_populations"]
    mothers = sources * m["mother_plants_sampled_per_population"]
    planned = mothers*m["enclosures"]*m["descendant_rosettes_per_mother_per_enclosure"]
    if (sources, mothers, planned) != (5, 50, 800):
        raise AssertionError("Original sample design no longer matches Methods")
    if m["inland_garden_is_original_seed_collection_site"] or not m[
        "inland_enclosure_near_inland_population_1_but_not_at_original_inland_sites"
    ]:
        raise AssertionError("Inland garden incorrectly promoted to original source site")
    return {
        "schema":"ttf_Arany2009_original_five_source_two_garden_maternal_grain_v01",
        "original_article_doi":DOI,
        "status":"FIVE_SOURCE_POPULATIONS_TWO_ENCLOSURES_NO_INDEPENDENT_GPS_OR_LGM_PANEL",
        "original_dune_seed_source_populations":m["dune_original_seed_populations"],
        "original_inland_seed_source_populations":m["inland_original_seed_populations"],
        "original_seed_source_population_labels":sources,
        "original_source_habitat_classes":m["distinct_original_source_habitat_classes"],
        "original_maternal_plants_sampled":mothers,
        "maternal_plants_per_population":m["mother_plants_sampled_per_population"],
        "experimental_garden_enclosures":m["enclosures"],
        "descendant_rosettes_per_mother_per_garden":m[
            "descendant_rosettes_per_mother_per_enclosure"],
        "nominal_transplanted_rosettes":planned,
        "nominal_rosettes_not_independent_original_source_GPS":True,
        "published_dune_pop1_pop2_distance_m":m[
            "dune_population_1_to_2_reported_distance_m"],
        "original_inland_seed_collections_occurred_at_experimental_inland_garden":False,
        "dune_garden_near_third_dune_seed_population":True,
        "original_source_population_geographic_precision_verified":0,
        "independent_climate_exposure_points_qualified":0,
        "climate_grid_assignment_allowed":False,
        "hargreaves_source_ID_labels":r["distinct_original_source_ID_labels"],
        "hargreaves_garden_ID_labels":r["distinct_original_garden_ID_labels"],
        "hargreaves_candidate_cells_not_independent_experiments":r[
            "source_only_candidate_same_garden_cells"],
        "original_methods_are_not_original_transplant_response_rows":True,
        "no_fitness_genetics_LGM_adaptation_plasticity_speciation_response_opened":True,
        "prior_Bontrager_44of50_stop_unchanged":True,
    }


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("--output",required=True,type=Path)
    args=parser.parse_args()
    d=json.loads(SHORTLIST.read_text(encoding="utf-8"))
    result=evaluate(d)
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(result,sort_keys=True,indent=2)+"\n",
                           encoding="utf-8")
    print(json.dumps(result,sort_keys=True))


if __name__ == "__main__":
    main()
