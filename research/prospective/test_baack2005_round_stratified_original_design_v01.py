"""Known-truth original Baack 2005 Table 1 design tests, outcome-blind."""
import copy,json,sys,unittest
from pathlib import Path
from unittest.mock import patch
sys.path.insert(0,str(Path(__file__).resolve().parent))
import baack2005_round_stratified_original_design_v01 as audit

class BaackStudySourceTests(unittest.TestCase):
    def setUp(self):self.d=json.loads(audit.SHORTLIST.read_text("utf-8"))
    def test_three_rounds_16_design_arms_each(self):
        a=audit.evaluate(self.d)
        self.assertEqual(len(a["seedling_experiment_rounds"]),3)
        self.assertEqual([x["directed_transplant_arms"] for x in a["seedling_experiment_rounds"]],[16,16,16])
    def test_six_sources_eight_gardens_and_ten_total_study_populations(self):
        a=audit.evaluate(self.d)
        self.assertEqual((a["original_seedling_source_populations"],
                          a["original_seedling_garden_populations"],
                          a["original_study_unique_populations_reported"]),(6,8,10))
    def test_48_round_arms_but_only_38_pairs(self):
        a=audit.evaluate(self.d)
        self.assertEqual((a["round_specific_directed_transplant_arms"],
                          a["distinct_source_garden_pairs_after_ignoring_round"]),(48,38))
        self.assertFalse(a["complete_6x8_factorial_present"])
    def test_10_home_arms_split_4_2_4(self):
        a=audit.evaluate(self.d)
        self.assertEqual([x["home_arms"] for x in a["seedling_experiment_rounds"]],[4,2,4])
        self.assertEqual(a["home_arms_with_round_preserved"],10)
    def test_cross_cytotype_environment_arms(self):
        self.assertEqual(audit.evaluate(self.d)["cross_cytotype_environment_arms_with_round_preserved"],24)
    def test_original_shortlist_cells_are_not_round_arms(self):
        a=audit.evaluate(self.d)
        self.assertEqual(a["source_only_hargreaves_candidate_cells"],8)
        self.assertTrue(a["Hargreaves_cells_not_directly_equivalent_to_round_stratified_arms"])
    def test_no_biological_response_and_no_GPS_promotion(self):
        a=audit.evaluate(self.d)
        self.assertEqual(a["source_GPS_originally_certified"],0)
        self.assertFalse(a["fitness_genetic_adaptation_plasticity_speciation_or_LGM_responses_opened"])
        self.assertTrue(a["old_bontrager_44_of_50_stop_unchanged"])
    def test_unqualified_data_cannot_silently_advance(self):
        d=copy.deepcopy(self.d)
        k=next(x for x in d["record_list"] if x["original_Hargreaves_study_label"]==audit.LABEL)
        k["historical_climate_contrast_qualified"]=True
        with self.assertRaisesRegex(ValueError,"qualification unexpectedly changed"):audit.evaluate(d)
    def test_original_table_round_drift_rejected(self):
        changed=copy.deepcopy(audit.ROUNDS)
        changed["Middle"]["garden"][0]="T-HPWM"
        with patch.dict(audit.ROUNDS,changed,clear=True):
            with self.assertRaises((ValueError,AssertionError)):audit.evaluate(self.d)

if __name__=="__main__":
    unittest.main()
