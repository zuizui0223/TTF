"""Frozen sensitivity MIROC source-only archive and response safeguards."""
import copy,json,sys,unittest
from pathlib import Path
from unittest.mock import patch
sys.path.insert(0,str(Path(__file__).resolve().parent))
import albano2025_MIROC_sensitivity_LGM_raster_unlabeled_v01 as a
class MIROCTest(unittest.TestCase):
 @classmethod
 def setUpClass(cls):cls.d=json.loads(a.CONTRACT.read_text("utf-8"))
 def test_original_miroc_secondary_role_not_primary(self):
  self.assertEqual((self.d["models"][1]["code"],self.d["models"][1]["role"]),("MIROC-ESM","sensitivity"))
 def test_wrong_source_link_rejected(self):
  d=copy.deepcopy(self.d);d["models"][1]["official_url"]="https://fake.test"
  with self.assertRaisesRegex(ValueError,"GCM role/source changed"):a.evaluate(d,inspector=lambda x:{})
 def test_3_failed_index_requests_do_not_open_any_biological_response(self):
  def fail(u):raise ValueError("network outage")
  x=a.evaluate(self.d,inspector=fail)
  self.assertEqual(len(x["source_events"]),3)
  self.assertFalse(x["LGM_values_extracted"])
  self.assertFalse(x["biological_response_read"])
 def test_valid_original_miroc_metadata_has_19_members(self):
  self.assertEqual(set(a.NAMES),{"BIO1","BIO12"})
  self.assertTrue(all(name.startswith("mrlgmbi") for name in a.NAMES.values()))
if __name__=="__main__":unittest.main()
