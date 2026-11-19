"""Original Melo 2014 coordinate-vs-prose source-only audit tests."""
import copy
import json
import sys
import unittest
from pathlib import Path
from unittest.mock import patch

sys.path.insert(0,str(Path(__file__).resolve().parent))
import melo2014_original_coordinate_distance_gate_v01 as audit

class MeloOriginalSourceTests(unittest.TestCase):
    def setUp(self):
        self.d = json.loads(audit.SHORTLIST.read_text(encoding="utf-8"))

    def test_original_dms_focal_pair(self):
        r = audit.evaluate(self.d)
        self.assertEqual(r["source_populations_for_focal_field_transplant"],2)
        self.assertTrue(r["publisher_coordinates_DMS_syntax_valid"])

    def test_great_circle_distance_known_value(self):
        r = audit.evaluate(self.d)
        self.assertAlmostEqual(r["great_circle_distance_km_spherical"],3.505057,places=5)
        self.assertAlmostEqual(r["distance_difference_km"],2.705057,places=5)

    def test_text_distance_discrepancy(self):
        r = audit.evaluate(self.d)
        self.assertEqual(r["text_approx_separation_km"],0.8)
        self.assertAlmostEqual(r["coordinate_to_text_distance_ratio"],4.381322,places=5)

    def test_no_gps_certification_after_published_coordinates(self):
        r = audit.evaluate(self.d)
        self.assertEqual(r["independently_verified_original_source_GPS"],0)
        self.assertFalse(r["historical_climate_raster_extraction_qualified"])
        self.assertFalse(r["potential_coordinate_error_proven"])

    def test_two_focal_sources_not_ten_independent_trials(self):
        r = audit.evaluate(self.d)
        self.assertEqual(r["hargreaves_candidate_comparison_cells_not_independent_experiments"],10)
        self.assertEqual(r["original_hargreaves_source_ID_labels"],2)
        self.assertFalse(r["single_paper_standalone_modern_plus_LGM_effect_estimable"])

    def test_no_fitness_genetics_or_historical_response(self):
        r = audit.evaluate(self.d)
        self.assertFalse(r["fitness_genetics_adaptation_plasticity_speciation_LGM_response_opened"])
        self.assertTrue(r["old_Bontrager_44_of_50_stop_unchanged"])

    def test_invalid_dms_seconds(self):
        with self.assertRaisesRegex(ValueError,"out of range"):
            audit.parse_dms("28°19′61.00″S","latitude")

    def test_hemisphere_mismatch(self):
        with self.assertRaisesRegex(ValueError,"hemisphere"):
            audit.parse_dms("153°34′46.82″S","longitude")

    def test_source_shortlist_gate_fails_closed(self):
        d = copy.deepcopy(self.d)
        item = next(x for x in d["record_list"] if x["original_Hargreaves_study_label"] == audit.LABEL)
        item["historical_climate_contrast_qualified"] = True
        with self.assertRaisesRegex(ValueError,"admission gate changed"):
            audit.evaluate(d)

    def test_publication_coordinate_drift_detected(self):
        bad = copy.deepcopy(audit.ORIGINAL_METHOD_COORDINATES)
        bad["Dune"]["latitude"] = "28°21′19.66″S"
        with patch.dict(audit.ORIGINAL_METHOD_COORDINATES,bad,clear=True):
            # Primary-vs-thesis agreement is checked before the distance ratio.
            with self.assertRaisesRegex(ValueError,"Cross-paper published coordinates no longer concordant"):
                audit.evaluate(self.d)

    def test_companion_thesis_separates_cabarita_and_hastings(self):
        r = audit.evaluate(self.d)
        x = r["published_site_name_and_seconds_concordance"]
        self.assertEqual(x["Dune"]["thesis_place_name"], "Cabarita Beach")
        self.assertEqual(x["Dune"]["thesis_population_code"], "SD03")
        self.assertEqual(x["Headland"]["thesis_place_name"], "Hastings Point")
        self.assertEqual(x["Headland"]["thesis_population_code"], "RH02")

    def test_crosspaper_DMS_within_one_arcsecond_not_same_cohort(self):
        r = audit.evaluate(self.d)
        x = r["published_site_name_and_seconds_concordance"]
        self.assertAlmostEqual(x["Dune"]["latitude_difference_arcsec"],0.66,places=5)
        self.assertAlmostEqual(x["Dune"]["longitude_difference_arcsec"],0.04,places=5)
        self.assertAlmostEqual(x["Headland"]["latitude_difference_arcsec"],0.07,places=5)
        self.assertAlmostEqual(x["Headland"]["longitude_difference_arcsec"],0.82,places=5)
        self.assertFalse(r["same_original_transplant_individuals_crosspaper_verified"])
        self.assertTrue(r["crosspaper_site_name_match_is_not_independent_original_GPS"])

    def test_companion_thesis_site_coordinate_drift_fails_closed(self):
        changed=copy.deepcopy(audit.THESIS_SITE_REFERENCE)
        changed["Headland"]["latitude"]="28°20′45″S"
        with patch.dict(audit.THESIS_SITE_REFERENCE,changed,clear=True):
            with self.assertRaisesRegex(ValueError,"no longer concordant"):
                audit.evaluate(self.d)

    def test_haversine_identity(self):
        self.assertEqual(audit.haversine((0,0),(0,0)),0)

if __name__=="__main__":
    unittest.main()
