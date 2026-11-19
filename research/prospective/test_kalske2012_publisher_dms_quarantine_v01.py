"""Kalske publisher DMS syntax and quarantine known-truth tests."""
import sys,unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parent))
from kalske2012_publisher_dms_quarantine_v01 import evaluate,validate
class DmsTests(unittest.TestCase):
    def test_one_syntax_valid(self):
        x=evaluate();self.assertEqual(x["valid_publisher_DMS_pairs"],1)
        self.assertEqual(x["certified_original_source_GPS"],0)
    def test_minute_out_of_range(self):
        self.assertIn("MINUTES_OUT_OF_RANGE",validate("15° 63′ 73″ E","longitude"))
    def test_second_out_of_range(self):
        self.assertIn("SECONDS_OUT_OF_RANGE",validate("15° 52′ 62″ E","longitude"))
    def test_valid_does_not_certify_source_gps(self):
        x=evaluate();u=next(z for z in x["entries"] if z["site"]=="Utö")
        self.assertTrue(u["valid_dms_pair"]);self.assertFalse(u["certified_source_GPS"])
    def test_no_cross_study_substitution(self):
        self.assertTrue(evaluate()["no_coordinate_correction_or_cross_study_substitution"])
    def test_no_biological_response(self):
        self.assertTrue(evaluate()["no_fitness_genetics_LGM_or_adaptation_response_opened"])
    def test_malformed_string_fails(self):
        self.assertEqual(validate("60 14.059 N","latitude"),["FORMAT_INVALID"])
if __name__=="__main__":unittest.main()
