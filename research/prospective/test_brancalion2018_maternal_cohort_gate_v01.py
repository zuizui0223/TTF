"""Original maternal cohort source-only regression tests for Brancalion 2018."""
import copy,json,sys,unittest
from pathlib import Path
from unittest.mock import patch
sys.path.insert(0,str(Path(__file__).resolve().parent))
import brancalion2018_maternal_cohort_gate_v01 as src

class BrancalionSourceTests(unittest.TestCase):
    def setUp(self):
        self.s=json.loads(src.SHORTLIST.read_text(encoding="utf-8"))

    def study(self,s=None):
        r=(s or self.s)["record_list"]
        return next(x for x in r if x["original_Hargreaves_study_label"]==src.LABEL)

    def test_three_populations_and_thirty_mothers(self):
        r=src.evaluate(self.s)
        self.assertEqual(r["published_original_source_populations"],3)
        self.assertEqual(r["published_mother_palms_per_population"],10)
        self.assertEqual(r["published_distinct_mother_palms_total"],30)

    def test_three_modalities_share_one_source_cohort(self):
        r=src.evaluate(self.s)
        self.assertEqual(r["experiment_modality_count"],3)
        self.assertTrue(r["modalities_share_original_maternal_seed_lots"])
        self.assertEqual(r["independent_maternal_source_cohorts_across_modalities"],1)

    def test_eight_candidate_cells_are_not_eight_trials(self):
        r=src.evaluate(self.s)
        self.assertEqual(r["hargreaves_same_garden_candidate_cells_not_independent"],8)
        self.assertEqual(r["hargreaves_original_source_ID_labels"],3)

    def test_no_new_qualified_panel_or_outcome(self):
        r=src.evaluate(self.s)
        self.assertEqual(r["source_population_GPS_independently_verified"],0)
        self.assertFalse(r["experiment_cohort_IDs_independently_certified"])
        self.assertFalse(r["read_fitness_SNP_LGM_or_adaptation_response"])
        self.assertTrue(r["bontrager_44_of_50_stop_unchanged"])

    def test_study_geometry_drift_rejected(self):
        s=copy.deepcopy(self.s)
        self.study(s)["distinct_original_source_ID_labels"]=4
        with self.assertRaisesRegex(ValueError,"source geometry changed"):
            src.evaluate(s)

    def test_duplicate_study_rejected(self):
        s=copy.deepcopy(self.s)
        s["record_list"].append(copy.deepcopy(self.study(s)))
        with self.assertRaisesRegex(ValueError,"not exactly one"):
            src.evaluate(s)

    def test_unreviewed_promotion_rejected(self):
        s=copy.deepcopy(self.s)
        self.study(s)["primary_experiment_identity_certified"]=True
        with self.assertRaisesRegex(ValueError,"admission unexpectedly changed"):
            src.evaluate(s)

    def test_maternal_math_drift_rejected(self):
        with patch.dict(src.PUBLISHED_METHODS,{"mother_palms_total":31}):
            with self.assertRaisesRegex(AssertionError,"maternal source count"):
                src.evaluate(self.s)

if __name__=="__main__":
    unittest.main()
