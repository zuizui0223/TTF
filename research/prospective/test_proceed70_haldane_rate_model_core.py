"""Synthetic-only Haldane rate model mathematics; no real PROCEED data."""
import sys,unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parent))
import numpy as np
from proceed70_haldane_rate_model_core import (
    ridge_predict,equal_species_weights,known_truth_labels,
    species_mse_comparison,wilson95,SEED
)

class RatesCoreTests(unittest.TestCase):
    def setUp(self):
        rng=np.random.default_rng(20261114)
        self.X=rng.normal(size=(120,21))
        self.H=rng.normal(size=120)
        self.sp=np.array([f"Species_{k//5}" for k in range(120)])
        self.st=np.array([f"Study_{k//5}" for k in range(120)])
        self.roles=np.array(["train" if k<60 else "heldout" for k in range(120)])
        self.rng=rng

    def test_weighted_ridge_design_and_single_history(self):
        y=self.rng.normal(size=120)
        p0=ridge_predict(self.X,y,self.sp,self.roles)
        p1=ridge_predict(self.X,y,self.sp,self.roles,history=self.H)
        self.assertEqual(p0.shape,(120,))
        self.assertTrue(np.isfinite(p1).all())
        with self.assertRaises(ValueError):
            ridge_predict(self.X,y,self.sp,self.roles,alpha=0.3)
        with self.assertRaises(ValueError):
            ridge_predict(self.X[:,:20],y,self.sp,self.roles)

    def test_training_only_fit_not_affected_by_holdout_outcomes(self):
        y=self.rng.normal(size=120)
        a=ridge_predict(self.X,y,self.sp,self.roles,history=self.H)
        y[60:]=1e5
        b=ridge_predict(self.X,y,self.sp,self.roles,history=self.H)
        np.testing.assert_allclose(a,b,rtol=0,atol=0)

    def test_original_species_equal_weights(self):
        w=equal_species_weights(self.sp,self.roles)
        self.assertAlmostEqual(w.sum(),1)
        self.assertEqual(float(w[60:].sum()),0)
        self.assertAlmostEqual(float(w[:5].sum()),float(w[5:10].sum()))

    def test_frozen_known_truth_and_species_bootstrap(self):
        y=known_truth_labels(self.X,self.H,self.sp,self.st,self.roles,
                             np.random.default_rng(SEED),0)
        a=ridge_predict(self.X,y,self.sp,self.roles)
        same=species_mse_comparison(y,a,a,self.sp,self.roles,reps=199,seed=SEED)
        self.assertEqual(same["relative_MSE_gain"],0)
        self.assertFalse(same["positive"])
        perfect=species_mse_comparison(y,a,y,self.sp,self.roles,reps=199,seed=SEED)
        self.assertGreater(perfect["relative_MSE_gain"],0.9)
        self.assertGreater(perfect["species_bootstrap_9833"][0],0.9)

    def test_source_simulations_have_finite_effect_and_wilson_guard(self):
        r=known_truth_labels(self.X,self.H,self.sp,self.st,self.roles,
                             np.random.default_rng(SEED),1)
        self.assertEqual(r.shape,(120,))
        self.assertTrue(np.isfinite(r).all())
        self.assertLess(wilson95(0,100)[1],0.1)
        self.assertGreater(wilson95(100,100)[0],0.8)

if __name__=="__main__":unittest.main()
