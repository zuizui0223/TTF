"""Fail-closed two-study scope and 10-paper admission closure regression tests."""
import copy,json,sys,unittest
from pathlib import Path
from unittest.mock import patch
sys.path.insert(0,str(Path(__file__).resolve().parent))
import volis_kindell_original_study_scope_closure_v01 as study

class FinalTwoShortlistTests(unittest.TestCase):
    def setUp(self):
        self.d=json.loads(study.SHORTLIST.read_text(encoding="utf-8"))

    def test_ten_papers_two_scoped(self):
        x=study.evaluate(self.d)
        self.assertEqual(x["original_shortlist_label_count"],10)
        self.assertEqual(x["newly_scoped_last_two_original_studies"],2)
        self.assertEqual(len(x["studies"]),2)

    def test_volis_two_source_sites_original_methods(self):
        x=study.evaluate(self.d)["studies"]
        volis=next(s for s in x if s["study_label"].startswith("Volis"))
        self.assertEqual(volis["original_source_populations"],2)
        self.assertEqual(volis["original_field_garden_sites"],2)
        self.assertIn("AUTHOR_HOSTED_ORIGINAL",volis["evidence_resolution"])

    def test_volis_two_maternal_pool_sizes_but_not_extra_populations(self):
        v=next(s for s in study.evaluate(self.d)["studies"] if s["study_label"].startswith("Volis"))
        d=v["published_source_only_design_traits"]
        self.assertEqual(d["cleared_vegetation_seed_pool_parents_per_origin"],3)
        self.assertEqual(d["intact_vegetation_seedling_parents_per_origin"],10)
        self.assertEqual(v["original_source_populations"],2)

    def test_kindell_two_sources_not_six_geographic_gardens(self):
        k=next(s for s in study.evaluate(self.d)["studies"] if s["study_label"].startswith("Kindell"))
        self.assertEqual((k["original_source_populations"],k["original_field_garden_sites"]),(2,2))
        self.assertEqual(k["Hargreaves_garden_labels_NOT_certified_as_distinct_field_geographies"],6)

    def test_kindell_only_published_abstract_verified(self):
        x=study.evaluate(self.d)
        self.assertFalse(x["full_original_methods_verified_for_kindell"])
        k=next(s for s in x["studies"] if s["study_label"].startswith("Kindell"))
        self.assertEqual(k["published_source_only_design_traits"]["original_age_treatments"],3)

    def test_no_gps_or_history_qualification(self):
        x=study.evaluate(self.d)
        self.assertEqual(x["existing_shortlist_verified_original_source_GPS_studies"],0)
        self.assertEqual(x["existing_shortlist_independent_original_experiments_certified"],0)
        self.assertEqual(x["historical_source_contrast_qualified"],0)
        self.assertFalse(x["independent_10_paper_panel_qualified"])

    def test_no_outcomes_or_frozen_science_modified(self):
        x=study.evaluate(self.d)
        self.assertTrue(x["no_fitness_genetics_LGM_adaptation_plasticity_speciation_response_opened"])
        self.assertTrue(x["prior_Bontrager_44of50_stop_unchanged"])

    def test_source_label_drift_rejected(self):
        d=copy.deepcopy(self.d)
        row=next(s for s in d["record_list"] if s["original_Hargreaves_study_label"]=="Volisetal_Oecologia_2002")
        row["distinct_original_source_ID_labels"]=4
        with self.assertRaisesRegex(ValueError,"grain mismatch"):study.evaluate(d)

    def test_garden_label_drift_rejected(self):
        d=copy.deepcopy(self.d)
        row=next(s for s in d["record_list"] if s["original_Hargreaves_study_label"]=="Kindelletal_JEcol_1996")
        row["distinct_original_garden_ID_labels"]=2
        with self.assertRaisesRegex(ValueError,"grain mismatch"):study.evaluate(d)

    def test_unverified_outcome_promotion_rejected(self):
        d=copy.deepcopy(self.d)
        row=next(s for s in d["record_list"] if s["original_Hargreaves_study_label"]=="Kindelletal_JEcol_1996")
        row["historical_climate_contrast_qualified"]=True
        with self.assertRaisesRegex(ValueError,"qualification changed"):study.evaluate(d)

if __name__=="__main__":
    unittest.main()
