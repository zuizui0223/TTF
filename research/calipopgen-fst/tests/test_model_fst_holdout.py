import sys,unittest,numpy as np
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from model_fst_holdout import evaluate
def world(beta,seed):
 rng=np.random.default_rng(seed);g=[]
 for i in range(100):
  n=12
  x=rng.normal(size=(n,5))
  y=.2*x[:,0]+.1*x[:,1]+beta*x[:,4]+rng.normal(0,.03,n)
  x-=x.mean(axis=0);y-=y.mean()
  g.append({"species":str(i),"role":i%2,"x":x,"y":y})
 return g
class Tests(unittest.TestCase):
 def test_large_known_effect(self):
  r=evaluate(world(.4,213),bootstrap_draws=101)
  self.assertGreater(r["relative_heldout_MSE_gain"],.5)
  self.assertTrue(r["positive_gain_95_lower_gt_zero"])
 def test_null_does_not_look_positive(self):
  r=evaluate(world(0.,987),bootstrap_draws=101)
  self.assertFalse(r["positive_gain_95_lower_gt_zero"])
 def test_species_leakage_stops(self):
  x=np.zeros((5,5));y=np.zeros(5)
  with self.assertRaises(ValueError):
   evaluate([{"species":"A","role":0,"x":x,"y":y},{"species":"A","role":1,"x":x,"y":y}],bootstrap_draws=5)
if __name__=="__main__":unittest.main()
