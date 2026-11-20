"""Offline current baseline retry and no climate-response tests."""
import json,sys,unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parent))
import worldclim_matched_baseline_index_retry_v01 as audit

class BaselineArchiveTests(unittest.TestCase):
 @classmethod
 def setUpClass(cls):cls.contract=json.loads(audit.CONTRACT.read_text("utf-8"))
 def test_official_current_baseline_url(self):
  r=audit.evaluate(self.contract,worker=lambda u:{"status":"OFFICIAL_ZIP_CENTRAL_DIRECTORY_VERIFIED"})
  self.assertEqual(r["status"],"OFFICIAL_CURRENT_BASELINE_CENTRAL_INDEX_VERIFIED")
  self.assertIn("bio_2-5m_bil.zip",r["source_url"])
 def test_error_attempts_failed_closed(self):
  def fail(url):raise ValueError("Synthetic source failure")
  x=audit.evaluate(self.contract,worker=fail)
  self.assertEqual(len(x["attempts"]),3)
  self.assertEqual(x["status"],"CURRENT_BASELINE_INDEX_ACCESS_HOLD")
  self.assertIsNone(x["index"])
 def test_source_only_boundary(self):
  r=audit.evaluate(self.contract,worker=lambda u:{"status":"OFFICIAL_ZIP_CENTRAL_DIRECTORY_VERIFIED"})
  self.assertFalse(r["LGM_or_modern_raster_values_read"])
  self.assertFalse(r["fitness_or_genetic_response_read"])
if __name__=="__main__":unittest.main()
