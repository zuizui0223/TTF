import unittest
import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import numpy as np
from spatial_null_core import (standardize_predictors,draw_synthetic_labels,species_auc,
    training_history_unique_fraction,evaluate_species_holdout)
class SyntheticTests(unittest.TestCase):
 def fixture(self,ns=24,n=90):
    rng=np.random.default_rng(515);s=np.repeat(np.arange(ns),n)
    r=s%2
    xy=np.column_stack([-110+rng.normal(size=ns*n),45+rng.normal(size=ns*n)])
    z=np.column_stack([rng.normal(size=(ns*n,4)),rng.normal(size=(ns*n,4)),
                        np.abs(rng.normal(size=(ns*n,5)))])
    return z,xy,s,r,{i:10 for i in range(ns)}
 def test_training_scaling_never_sees_heldout(self):
    z,xy,s,r,c=self.fixture()
    a,mu,sd,clim=standardize_predictors(z,xy,s,r)
    z[r==1]+=1e4
    b,mu2,sd2,clim2=standardize_predictors(z,xy,s,r)
    np.testing.assert_array_equal(mu,mu2)
    np.testing.assert_array_equal(sd,sd2)
    np.testing.assert_array_equal(clim,clim2)
 def test_synthetic_cases_preserve_locked_count(self):
    z,xy,s,r,c=self.fixture()
    p,*_=standardize_predictors(z,xy,s,r)
    y=draw_synthetic_labels(p,xy,s,c,np.random.default_rng(11),.9)
    self.assertEqual(y.sum(),sum(c.values()))
    self.assertTrue(all(y[s==k].sum()==v for k,v in c.items()))
 def test_rank_auc_respects_ties(self):
    y=[1,0,1,0,1,0];v=[.9,.1,.7,.2,.2,.7]
    self.assertAlmostEqual(species_auc(y,v),7/9,places=8)
 def test_historical_redundancy_gate(self):
    z,xy,s,r,c=self.fixture()
    p,*_=standardize_predictors(z,xy,s,r)
    self.assertGreater(training_history_unique_fraction(p,s,r),.05)
    p[:,11]=p[:,0]+2*p[:,1]
    self.assertLess(training_history_unique_fraction(p,s,r),1e-9)
 def test_known_positive_and_species_disjoint_evaluation(self):
    z,xy,s,r,c=self.fixture()
    p,*_=standardize_predictors(z,xy,s,r)
    y=draw_synthetic_labels(p,xy,s,c,np.random.default_rng(22),3.0,0.)
    v=evaluate_species_holdout(p,y,s,r,99,77)
    self.assertEqual(v["train_species"],12)
    self.assertEqual(v["heldout_species"],12)
    self.assertGreater(v["mean_delta_auc"],0.)
    self.assertTrue(np.isfinite(v["ci_9833"]).all())
if __name__=="__main__":unittest.main()
