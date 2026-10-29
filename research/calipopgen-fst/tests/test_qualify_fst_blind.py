import sys,unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from qualify_fst_blind import frozen_binomial,gc_km
class TestGeoQual(unittest.TestCase):
 def test_binomial(self):
  self.assertEqual(frozen_binomial("Puma concolor californica"),"puma concolor")
 def test_haversine(self):
  self.assertAlmostEqual(gc_km((0,0),(0,1)),111.195,places=2)
 def test_bad_taxon(self):
  self.assertIsNone(frozen_binomial("Puma sp."))
if __name__=="__main__":unittest.main()
