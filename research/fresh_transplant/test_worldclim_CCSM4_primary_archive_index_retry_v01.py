"""Frozen primary CCSM4 archive identity and source-only retry tests."""
import json,sys,unittest,copy
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parent))
import worldclim_CCSM4_primary_archive_index_retry_v01 as a
class CCSM4Tests(unittest.TestCase):
 @classmethod
 def setUpClass(cls):cls.d=json.loads(a.CONTRACT.read_text("utf-8"))
 def test_precommitted_primary_role(self):
  x=a.evaluate(self.d,inspector=lambda u:{"status":"OFFICIAL_ZIP_CENTRAL_DIRECTORY_VERIFIED"})
  self.assertEqual(x["status"],"PRIMARY_CCSM4_ZIP_CENTRAL_INDEX_VERIFIED")
  self.assertIn("/cclgmbi_2-5m.zip",x["source_url"])
 def test_frozen_primary_cannot_substitute_MIROC(self):
  d=copy.deepcopy(self.d);d["models"][0]["code"]="MIROC-ESM"
  with self.assertRaisesRegex(ValueError,"primary GCM role changed"):a.evaluate(d,inspector=lambda u:{})
 def test_transport_failure_never_promotes_GCM(self):
  def fail(url):raise ValueError("Blocked")
  x=a.evaluate(self.d,inspector=fail)
  self.assertEqual(len(x["attempts"]),3)
  self.assertFalse(x["original_CCSM4_raster_values_extracted"])
  self.assertFalse(x["fitness_or_genetic_responses_opened"])
  self.assertEqual(x["status"],"PRIMARY_CCSM4_ZIP_INDEX_NOT_YET_ACCESSIBLE")
if __name__=="__main__":unittest.main()
