import sys,unittest,numpy as np
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from climate_feature_gate import partial_fraction
class TestInfoGate(unittest.TestCase):
 def test_historical_noise_survives(self):
  rng=np.random.default_rng(33)
  groups=[]
  for i in range(10):
   x=rng.normal(size=(20,5));x-=x.mean(axis=0);groups.append(x)
  self.assertGreater(partial_fraction(groups),.6)
 def test_exact_redundancy_removed(self):
  rng=np.random.default_rng(34);groups=[]
  for i in range(10):
   x=rng.normal(size=(20,4));h=x[:,0]+2*x[:,1];v=np.column_stack([x,h]);v-=v.mean(axis=0);groups.append(v)
  self.assertLess(partial_fraction(groups),1e-12)
if __name__=="__main__":unittest.main()
