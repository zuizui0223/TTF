"""Published methods identity audit: likely shared Oct-2004 adult transplant.

Only author-reported experimental designs and previously frozen nonresponse
location receipts are inspected, not plant performance or source response.
"""
import json
import unittest
from pathlib import Path

HERE=Path(__file__).resolve().parent

class Shared2004AdultTransplantTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.x=json.loads((HERE/"raabova2007_pankova2014_probable_shared_adult_cohort_v01.json").read_text())
        cls.geo=json.loads((HERE/"source_registry/raabova2007_pankova2014_six_site_exploratory_v01.json").read_text())
    def test_original_methods_match_not_only_reused_population_names(self):
        a=self.x["target_2007"]["adult_transplant_methods"]
        b=self.x["followup_2014"]["adult_transplant_methods"]
        fields=("seeds_collected","seed_sown","greenhouse","seedlings","field_transplant",
                "population_sites","origins_per_site","juvenile_replicates_per_origin_site",
                "blocks_rows_per_site","interplant_distance_cm")
        for key in fields:
            self.assertEqual(a[key],b[key],key)
        self.assertEqual(a["field_transplant"],"October 2004")
        self.assertEqual(a["juvenile_replicates_per_origin_site"],25)
        self.assertEqual(a["origins_per_site"],6)
        self.assertEqual(a["population_sites"],6)
    def test_shared_cohort_probable_individual_identity_not_certified(self):
        r=self.x["methods_identity_assessment"]
        self.assertTrue(r["supports_probable_shared_2004_adult_cohort"])
        self.assertFalse(r["individual_label_ID_crosspaper_directly_verified"])
        self.assertFalse(r["shared_plants_confirmed_with_persistent_individual_tag_match"])
        self.assertFalse(r["independent_replication_certified"])
        self.assertTrue(self.x["target_2007"]["study_2007_other_components"].startswith("Two seed sowing"))
    def test_six_site_original_geo_receipt_immutable(self):
        geo=self.geo
        self.assertEqual(geo["external_later_paper_doi"],self.x["followup_2014"]["doi"])
        self.assertEqual(geo["target_primary_paper_doi"],self.x["target_2007"]["doi"])
        self.assertEqual([z["distance_km"] for z in geo["site_details"]],
                         [0.799619,1.812728,1.928084,0.734489,1.311237,1.92598])
        self.assertEqual(geo["descriptive_sites_within_1km"],2)
        self.assertTrue(geo["no_adaptation_plasticity_speciation_genetic_history_climate_response_opened"])
    def test_new_inference_cannot_unfreeze_old_science(self):
        x=self.x["boundaries"]
        self.assertTrue(x["no_new_independent_experiment_added"])
        self.assertTrue(x["bontrager_previous_44_of_50_frozen_stop_unchanged"])
        self.assertTrue(x["PR155_remains_draft"])
        self.assertTrue(self.geo["old_invalid_Sudova_2010_vs_2007_receipt_remains_invalid"])

if __name__=="__main__":
    unittest.main()
