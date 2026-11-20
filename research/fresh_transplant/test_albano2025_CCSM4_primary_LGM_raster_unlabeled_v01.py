"""Frozen primary CCSM4 source-only original archive and wrapper tests."""
import copy,json,sys,unittest
from pathlib import Path
from unittest.mock import patch
sys.path.insert(0,str(Path(__file__).resolve().parent))
import albano2025_CCSM4_primary_LGM_raster_unlabeled_v01 as a
class PrimaryCCSM4Tests(unittest.TestCase):
 @classmethod
 def setUpClass(cls):
  cls.manifest=json.loads(a.MANIFEST.read_text("utf-8"))
  cls.contract=json.loads(a.CONTRACT.read_text("utf-8"))
  cls.grid=json.loads(a.engine.CELLS.read_text("utf-8"))
 def test_exact_two_primary_lgm_members(self):
  self.assertEqual(len(self.manifest["members"]),2)
  self.assertEqual({v["filename"] for v in self.manifest["members"]},
                   {"cclgmbi1.tif","cclgmbi12.tif"})
 def test_full_primary_archive_identity_from_original_index(self):
  self.assertEqual(self.manifest["original_archive_total_bytes"],219610901)
  self.assertEqual(self.manifest["original_archive_member_count"],19)
  self.assertEqual(self.manifest["original_source_run"],44755277981)
 def test_prespecified_primary_role_guard(self):
  d=copy.deepcopy(self.contract)
  d["models"][0]["role"]="sensitivity"
  with self.assertRaisesRegex(ValueError,"precommitted CCSM4"):
   a.evaluate(self.manifest,d,grid=self.grid)
 def test_wrong_primary_raster_member_guard(self):
  d=copy.deepcopy(self.manifest)
  d["members"][0]["filename"]="cclgmbi5.tif"
  with self.assertRaisesRegex(ValueError,"original BIO variables"):
   a.evaluate(d,self.contract,grid=self.grid)
 def test_original_engine_values_and_globals_restored(self):
  old=(a.engine.URL,a.engine.ARCHIVE_BYTES,a.engine.ENTRIES)
  def fake(grid):
   self.assertEqual(a.engine.ARCHIVE_BYTES,219610901)
   self.assertEqual(a.engine.ENTRIES["BIO1"]["crc32"],"51b70e63")
   self.assertEqual(a.engine.ENTRIES["BIO12"]["crc32"],"7d838cc2")
   self.assertEqual(len(grid["cells"]),95)
   return {"source_population_count":96,"unlabeled_grid_cell_count":95,
           "fitness_genetics_adaptation_responses_opened":False,
           "LGM_independent_predictor_gain_tested":False}
  with patch.object(a.engine,"evaluate",side_effect=fake):
   out=a.evaluate(self.manifest,self.contract,grid=self.grid)
  self.assertEqual((a.engine.URL,a.engine.ARCHIVE_BYTES,a.engine.ENTRIES),old)
  self.assertEqual(out["GCM_role"],"primary_precommitted")
  self.assertFalse(out["LGM_independent_predictor_gain_tested"])
 def test_missing_source_url_guard(self):
  d=copy.deepcopy(self.manifest)
  d["original_archive_url"]="https://example.org/cc.zip"
  with self.assertRaisesRegex(ValueError,"source URL changed"):
   a.evaluate(d,self.contract,grid=self.grid)
if __name__=="__main__":unittest.main()
