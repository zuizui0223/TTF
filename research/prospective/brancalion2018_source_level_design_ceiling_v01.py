"""Brancalion 2018 original Methods: independent site-level climate design ceiling.

Only published source/garden experiment structure, not outcomes or GPS.
Counts are nominal planned placements and NEVER biological observations.
Every provenance is one of three sampled forest populations.
"""
from __future__ import annotations
import argparse,json
from pathlib import Path

from brancalion2018_maternal_cohort_gate_v01 import (
    DOI, LABEL, SHORTLIST, evaluate as maternal_gate
)

# DOI 10.1002/ece3.4248, original Methods sections 2.1, 2.4-2.7.
SOURCE_POPULATIONS=3
MOTHERS_PER_SOURCE=10
GARDENS=3
BLOCKS_PER_GARDEN=8
PLOTS_PER_BLOCK=3
SEEDS_PER_PROVENANCE_PLOT=30
SEEDLINGS_PER_PROVENANCE_PLOT=10
COMMON_GARDEN_SAPLINGS_PER_MOTHER=2

def evaluate(source):
    gate=maternal_gate(source)
    if gate["published_original_source_populations"] != SOURCE_POPULATIONS:
        raise ValueError("source-population count drift")
    if gate["published_mother_palms_per_population"] != MOTHERS_PER_SOURCE:
        raise ValueError("maternal provenance count drift")
    if not gate["modalities_share_original_maternal_seed_lots"]:
        raise ValueError("independent cohort gate drift")
    total_blocks=GARDENS*BLOCKS_PER_GARDEN
    total_plots=total_blocks*PLOTS_PER_BLOCK
    planted_seeds=total_plots*SEEDS_PER_PROVENANCE_PLOT
    planted_seedlings=total_plots*SOURCE_POPULATIONS*SEEDLINGS_PER_PROVENANCE_PLOT
    common_garden_saplings=SOURCE_POPULATIONS*MOTHERS_PER_SOURCE*COMMON_GARDEN_SAPLINGS_PER_MOTHER
    if (total_blocks,total_plots,planted_seeds,planted_seedlings,common_garden_saplings)!=(24,72,2160,2160,60):
        raise AssertionError("published method placement arithmetic drift")
    controls={
      "intercept_and_one_modern_variable_and_one_LGM": 3,
      "intercept_and_two_modern_variables_and_one_LGM": 4,
      "intercept_and_seven_controls_and_one_LGM": 9,
    }
    return {
       "schema":"ttf_Brancalion2018_original_source_level_design_ceiling_v01",
       "primary_paper_doi":DOI,
       "status":"THREE_ORIGINAL_SOURCE_SITES_CANNOT_STANDALONE_IDENTIFY_MULTIVARIATE_LGM_EFFECT",
       "original_source_population_geo_exposure_units":SOURCE_POPULATIONS,
       "original_garden_sites":GARDENS,
       "original_maternal_trees":SOURCE_POPULATIONS*MOTHERS_PER_SOURCE,
       "seed_and_seedling_modalities_share_original_maternal_seed_lots":True,
       "original_reciprocal_garden_blocks_nominal":total_blocks,
       "original_reciprocal_plots_nominal":total_plots,
       "seed_transplant_planted_seed_positions_nominal":planted_seeds,
       "seedling_transplant_planted_seedling_positions_nominal":planted_seedlings,
       "common_garden_sapling_positions_nominal":common_garden_saplings,
       "all_placements_are_DESIGN_COUNTS_not_verified_observations":True,
       "max_source_level_climate_design_matrix_rank":SOURCE_POPULATIONS,
       "example_covariate_columns":controls,
       "intercept_plus_1modern_plus_LGM_would_have_residual_df_at_most":0,
       "intercept_plus_2modern_plus_LGM_must_be_rank_deficient":True,
       "intercept_plus_7controls_plus_LGM_must_be_rank_deficient":True,
       "multi_study_hierarchical_pooling_not_ruled_out":True,
       "published_forest_centroids_cannot_replace_original_mother_GPS":True,
       "independent_original_population_GPS_certified":0,
       "no_fitness_genetics_LGM_adaptation_or_speciation_response_opened":True,
       "bontrager_44_of_50_frozen_gate_unchanged":True,
       "interpretation":"Many offspring/blocks do not create independent source-level LGM exposures; no standalone residual degrees of freedom for intercept+modern+LGM on three sources.",
    }

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--output",type=Path,required=True)
    a=p.parse_args()
    src=json.loads(SHORTLIST.read_text(encoding="utf-8"))
    r=evaluate(src)
    a.output.parent.mkdir(parents=True,exist_ok=True)
    a.output.write_text(json.dumps(r,indent=2,sort_keys=True)+"\n",encoding="utf-8")
    print(json.dumps(r,sort_keys=True))
if __name__=="__main__":
    main()
