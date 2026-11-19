"""Regression tests for Kalske 2012 source-only DMS + geography v02 quarantine."""
import sys
import unittest
from pathlib import Path
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parent))
import kalske2012_dms_geographic_quarantine_v02 as audit


class KalskeGeographicTests(unittest.TestCase):
    def test_three_invalid_and_one_valid_dms(self):
        result = audit.evaluate()
        self.assertEqual((result["syntax_invalid_count"],
                          result["syntax_valid_but_geographically_inconsistent_count"]), (3, 1))

    def test_uto_syntax_valid_but_outside_region(self):
        uto = next(x for x in audit.evaluate()["sites"] if x["site"] == "Utö")
        self.assertTrue(uto["syntax_valid"])
        self.assertFalse(uto["inside_conservative_finland_envelope"])
        self.assertAlmostEqual(uto["converted_latitude"], 66.498333, places=6)
        self.assertAlmostEqual(uto["converted_longitude"], 15.357778, places=6)

    def test_zero_source_gps_certified_or_climate_extracted(self):
        result = audit.evaluate()
        self.assertEqual(result["certified_original_source_GPS_count"], 0)
        self.assertFalse(result["climate_grid_extraction_allowed"])
        self.assertFalse(result["fitness_genetics_adaptation_plasticity_speciation_LGM_response_opened"])

    def test_out_of_range_minute_fails_not_normalizes(self):
        val, err = audit.parse_strict_dms("15° 63′ 73″ E", "longitude")
        self.assertIsNone(val)
        self.assertIn("MINUTES_OUT_OF_RANGE", err)
        self.assertIn("SECONDS_OUT_OF_RANGE", err)

    def test_one_good_synthetic_coordinate_not_source_gps(self):
        site = audit.assess_site("synthetic", "60° 20′ 00″ N", "22° 00′ 00″ E")
        self.assertTrue(site["syntax_valid"])
        self.assertTrue(site["inside_conservative_finland_envelope"])
        self.assertFalse(site["original_source_GPS_certified"])
        self.assertEqual(site["admission"], "QUARANTINE")

    def test_hemisphere_mismatch_fails_closed(self):
        val, issues = audit.parse_strict_dms("15° 21′ 28″ N", "longitude")
        self.assertIsNone(val)
        self.assertIn("HEMISPHERE_MISMATCH", issues)

    def test_geographic_bound_strictness(self):
        for val, axis in (("90° 01′ 00″ N", "latitude"),
                          ("180° 00′ 01″ E", "longitude")):
            with self.subTest(val=val):
                _, issues = audit.parse_strict_dms(val, axis)
                self.assertIn("DEGREES_OUT_OF_RANGE", issues)

    def test_publisher_string_drift_rejected(self):
        with patch.dict(audit.PUBLISHED, {"Utö": ("60° 20′ 00″ N", "22° 00′ 00″ E")}):
            with self.assertRaisesRegex(ValueError, "re-audit needed"):
                audit.evaluate()


if __name__ == "__main__":
    unittest.main()
