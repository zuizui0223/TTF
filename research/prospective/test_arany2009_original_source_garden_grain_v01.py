"""Source-only Arany 2009 maternal-line versus site qualification tests."""
import copy
import json
import sys
import unittest
from pathlib import Path
from unittest.mock import patch

sys.path.insert(0,str(Path(__file__).resolve().parent))
import arany2009_original_source_garden_grain_v01 as study

class AranySourceTests(unittest.TestCase):
    def setUp(self):
        self.d=json.loads(study.SHORTLIST.read_text(encoding="utf-8"))

    def test_five_original_source_populations(self):
        a=study.evaluate(self.d)
        self.assertEqual((a["original_dune_seed_source_populations"],
                          a["original_inland_seed_source_populations"]), (3,2))
        self.assertEqual(a["original_seed_source_population_labels"],5)

    def test_50_source_mothers_but_no_extra_gps(self):
        a=study.evaluate(self.d)
        self.assertEqual((a["original_maternal_plants_sampled"],
                          a["maternal_plants_per_population"]),(50,10))
        self.assertEqual(a["original_source_population_geographic_precision_verified"],0)

    def test_800_rosettes_not_800_climate_exposures(self):
        a=study.evaluate(self.d)
        self.assertEqual(a["nominal_transplanted_rosettes"],800)
        self.assertEqual(a["experimental_garden_enclosures"],2)
        self.assertTrue(a["nominal_rosettes_not_independent_original_source_GPS"])

    def test_inland_garden_near_but_not_at_original_sources(self):
        a=study.evaluate(self.d)
        self.assertFalse(a["original_inland_seed_collections_occurred_at_experimental_inland_garden"])
        self.assertTrue(a["dune_garden_near_third_dune_seed_population"])

    def test_two_dune_origins_only_five_meters_apart(self):
        a=study.evaluate(self.d)
        self.assertEqual(a["published_dune_pop1_pop2_distance_m"],5)
        self.assertEqual(a["original_source_habitat_classes"],2)

    def test_secondary_six_candidate_cells_not_six_studies(self):
        a=study.evaluate(self.d)
        self.assertEqual(a["hargreaves_source_ID_labels"],5)
        self.assertEqual(a["hargreaves_garden_ID_labels"],2)
        self.assertEqual(a["hargreaves_candidate_cells_not_independent_experiments"],6)

    def test_no_outcome_or_climate_opening(self):
        a=study.evaluate(self.d)
        self.assertTrue(a["no_fitness_genetics_LGM_adaptation_plasticity_speciation_response_opened"])
        self.assertFalse(a["climate_grid_assignment_allowed"])
        self.assertTrue(a["prior_Bontrager_44of50_stop_unchanged"])

    def test_shortlist_promotion_fails_closed(self):
        d=copy.deepcopy(self.d)
        row=next(x for x in d["record_list"] if x["original_Hargreaves_study_label"]==study.LABEL)
        row["primary_experiment_identity_certified"]=True
        with self.assertRaisesRegex(ValueError,"qualification changed"):
            study.evaluate(d)

    def test_proxy_inland_garden_mislabel_fails_closed(self):
        with patch.dict(study.METHODS,{"inland_garden_is_original_seed_collection_site":True}):
            with self.assertRaisesRegex(AssertionError,"promoted"):
                study.evaluate(self.d)

if __name__=="__main__":
    unittest.main()
