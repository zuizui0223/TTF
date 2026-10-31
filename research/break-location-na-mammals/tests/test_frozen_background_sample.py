import sys,unittest,random
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from frozen_background_sample import select_background,hash_key
class Tests(unittest.TestCase):
 def test_independent_of_input_order(self):
  p=[(-110+i*.05,40+i*.00001) for i in range(1400)]
  a=select_background("Sciurus aberti",p)
  random.Random(177).shuffle(p)
  b=select_background("Sciurus aberti",p)
  self.assertEqual(a,b)
  self.assertEqual(len(a),1000)
 def test_keep_all_below_cap(self):
  self.assertEqual(len(select_background("Marmota broweri",[(-150.0+i*.01,60.) for i in range(33)])),33)
 def test_species_identity_is_used(self):
  p=[(-110+i*.05,40.) for i in range(1100)]
  self.assertNotEqual(select_background("Marmota broweri",p),select_background("Sciurus aberti",p))
 def test_disallows_posthoc_cap_change(self):
  with self.assertRaises(ValueError):
   select_background("Sciurus aberti",[(-110,40)],max_sites=500)
if __name__=="__main__":unittest.main()
