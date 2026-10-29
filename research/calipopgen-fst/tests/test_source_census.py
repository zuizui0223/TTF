import sys,unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from source_census import choose_article,normalize_detail

class TestCensus(unittest.TestCase):
 def test_exact_article(self):
  items=[{"id":123,"title":"CaliPopGen Dataset 1: Population Genetic Diversity"},
         {"id":456,"title":"CaliPopGen Dataset 2: Pairwise Population Differentiation"}]
  self.assertEqual(choose_article(items),456)
 def test_ambiguous_fails(self):
  with self.assertRaises(ValueError):choose_article([])
 def test_metadata_only(self):
  item={"id":456,"title":"CaliPopGen Dataset 2: Pairwise Population Differentiation",
       "doi":"10.x/y","files":[{"id":2,"name":"pairwise.csv","size":100,
                              "supplied_md5":"abc","GenDist":.98}]}
  x=normalize_detail(item)
  self.assertFalse(x["Fst_numeric_opened"])
  self.assertNotIn("0.98",repr(x))
if __name__=="__main__":unittest.main()
