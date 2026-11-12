"""Stage24 numerical and authority tests use only synthetic predictors and source rules."""
import sys
import unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))

import numpy as np
from spatial_null_core import standardize_predictors
from stage24_one_shot_break_auc import (
    frozen_original_predictor_scale_for_jitter,mandatory_authorities
)

class Stage24Tests(unittest.TestCase):
    def test_frozen_authority_is_source_only_and_complete(self):
        import json
        p=Path(__file__).resolve().parents[1]/"stage24_empirical_break_location_one_shot_contract_v01.json"
        c=json.loads(p.read_text())
        sc=mandatory_authorities(c)
        self.assertEqual(sc["bootstrap"]["draws_for_future_one_shot"],1999)
        self.assertTrue(sc["synthetic_gate"]["not_recalibrate_effect_size_or_cutoffs_after_results"])

    def test_original_training_scales_reproduce_identical_features(self):
        rng=np.random.default_rng(20261029)
        x=np.column_stack([rng.normal(size=(48,8))*4+14,
                           rng.uniform(0,200,size=(48,4)),
                           rng.uniform(0,30,size=48)])
        xy=np.column_stack([rng.uniform(-100,-90,48),rng.uniform(20,50,48)])
        sp=np.repeat(np.arange(4),12)
        role=np.repeat([0,0,1,1],12)
        original,mu,scale,sd=standardize_predictors(x,xy,sp,role)
        copied=frozen_original_predictor_scale_for_jitter(x.copy(),xy.copy(),mu,scale,sd)
        np.testing.assert_allclose(copied,original,rtol=0,atol=0)
        shifted=x.copy()
        shifted[:,0]+=2
        other=frozen_original_predictor_scale_for_jitter(shifted,xy,mu,scale,sd)
        self.assertTrue(np.isfinite(other).all())
        self.assertGreater(float(np.max(np.abs(other[:,-1]-original[:,-1]))),0)

    def test_source_contract_forbids_retuning(self):
        import json
        c=json.loads((Path(__file__).resolve().parents[1]/"stage24_empirical_break_location_one_shot_contract_v01.json").read_text())
        self.assertTrue(c["scope"]["no_subgroup_rescue"])
        self.assertTrue(c["scope"]["no_model_or_predictor_fitting_until_every_source_and_replay_SHA_gate_passes"])
        self.assertEqual(c["model"]["train_ids"],44)
        self.assertEqual(c["model"]["evaluation_ids"],36)

if __name__=="__main__":unittest.main()
