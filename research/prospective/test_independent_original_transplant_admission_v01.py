"""Fail-closed original-transplant source admission regression tests."""
import copy,json,sys,unittest
from pathlib import Path
HERE=Path(__file__).resolve().parent;sys.path.insert(0,str(HERE))
from independent_original_transplant_admission_v01 import SOURCE,evaluate
class AdmissionTests(unittest.TestCase):
    def setUp(self):self.d=json.loads(SOURCE.read_text())
    def test_ten_sources_zero_certified(self):
        r=evaluate(self.d)
        self.assertEqual(r["bibliographic_labels_checked"],10)
        self.assertEqual(r["candidate_cells_not_independent_replicates"],106)
        self.assertEqual(r["certified_independent_original_experiments"],0)
        self.assertEqual(r["certified_original_source_population_GPS_studies"],0)
    def test_no_biological_opening(self):
        r=evaluate(self.d)
        self.assertFalse(r["fitness_genetics_plasticity_speciation_or_LGM_response_opened"])
        self.assertTrue(r["prior_bontrager_44_of_50_stop_unchanged"])
    def test_missing_study_fails(self):
        d=copy.deepcopy(self.d);d["record_list"].pop()
        with self.assertRaisesRegex(ValueError,"set drift"):evaluate(d)
    def test_uncertified_study_not_auto_promoted(self):
        d=copy.deepcopy(self.d);d["record_list"][0]["primary_experiment_identity_certified"]=True
        with self.assertRaisesRegex(ValueError,"boundary changed"):evaluate(d)
    def test_new_gps_count_requires_reaudit(self):
        d=copy.deepcopy(self.d);d["summary"]["study_labels_with_original_source_lat_lon_verified"]=1
        with self.assertRaisesRegex(ValueError,"re-audit"):evaluate(d)
    def test_no_fitness_opened(self):
        d=copy.deepcopy(self.d);d["record_list"][1]["genetic_adaptation_plasticity_speciation_effect_opened"]=True
        with self.assertRaisesRegex(ValueError,"boundary changed"):evaluate(d)
    def test_no_study_cell_inflation(self):
        d=copy.deepcopy(self.d);d["record_list"][1]["source_only_candidate_same_garden_cells"]=-1
        with self.assertRaisesRegex(ValueError,"negative"):evaluate(d)
if __name__=="__main__":unittest.main()
