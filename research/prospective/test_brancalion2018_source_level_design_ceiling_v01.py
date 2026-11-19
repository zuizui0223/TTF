"""Known-truth tests for Brancalion 2018 climate-exposure rank ceiling."""
import copy,json,sys,unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parent))
from brancalion2018_source_level_design_ceiling_v01 import SHORTLIST,evaluate

class BrancalionDesignCeilingTests(unittest.TestCase):
    def setUp(self):self.source=json.loads(SHORTLIST.read_text("utf-8"))
    def test_exact_three_source_exposures(self):
        r=evaluate(self.source)
        self.assertEqual(r["original_source_population_geo_exposure_units"],3)
        self.assertEqual(r["max_source_level_climate_design_matrix_rank"],3)

    def test_mothers_not_sites(self):
        r=evaluate(self.source)
        self.assertEqual(r["original_maternal_trees"],30)
        self.assertNotEqual(r["original_maternal_trees"],r["original_source_population_geo_exposure_units"])

    def test_source_garden_block_and_plot_design(self):
        r=evaluate(self.source)
        self.assertEqual((r["original_garden_sites"],r["original_reciprocal_garden_blocks_nominal"],
                          r["original_reciprocal_plots_nominal"]),(3,24,72))

    def test_seed_positions_nominal(self):
        r=evaluate(self.source)
        self.assertEqual(r["seed_transplant_planted_seed_positions_nominal"],2160)
        self.assertTrue(r["all_placements_are_DESIGN_COUNTS_not_verified_observations"])

    def test_seedling_and_common_garden_positions_nominal(self):
        r=evaluate(self.source)
        self.assertEqual(r["seedling_transplant_planted_seedling_positions_nominal"],2160)
        self.assertEqual(r["common_garden_sapling_positions_nominal"],60)

    def test_no_residual_df_for_one_modern_plus_history(self):
        r=evaluate(self.source)
        self.assertEqual(r["example_covariate_columns"]["intercept_and_one_modern_variable_and_one_LGM"],3)
        self.assertEqual(r["intercept_plus_1modern_plus_LGM_would_have_residual_df_at_most"],0)

    def test_two_modern_variables_cannot_standalone_fit(self):
        r=evaluate(self.source)
        self.assertEqual(r["example_covariate_columns"]["intercept_and_two_modern_variables_and_one_LGM"],4)
        self.assertTrue(r["intercept_plus_2modern_plus_LGM_must_be_rank_deficient"])

    def test_no_outcomes_and_prior_stop_unmodified(self):
        r=evaluate(self.source)
        self.assertTrue(r["no_fitness_genetics_LGM_adaptation_or_speciation_response_opened"])
        self.assertTrue(r["bontrager_44_of_50_frozen_gate_unchanged"])
        self.assertEqual(r["independent_original_population_GPS_certified"],0)

    def test_original_source_gate_must_fail_closed(self):
        d=copy.deepcopy(self.source)
        item=next(x for x in d["record_list"]
                  if x["original_Hargreaves_study_label"]=="Brancalionetal_EcolEvol_2018")
        item["historical_climate_contrast_qualified"]=True
        with self.assertRaisesRegex(ValueError,"admission unexpectedly changed"):
            evaluate(d)

if __name__=="__main__":unittest.main()
