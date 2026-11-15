"""Synthetic-only Haldane 21-climate/method source design qualification tests."""
import sys,unittest,json
from pathlib import Path
from collections import Counter
sys.path.insert(0,str(Path(__file__).resolve().parent))
import numpy as np
from proceed70_haldane_source_null_calibration import (
    method_features,source_history_concentration,known_truth_gate,CONTRACT
)


class SyntheticSourceCalibrationTests(unittest.TestCase):
    def test_fixed_method_feature_sequence_and_lexical_source_categories(self):
        x=[("Plantus alpha","study1",(100.,30.),(101.,32.))]
        meta={x[0]:{"generations":[1.,2.,3.],
                    "trait":Counter({"size":2,"physio":1}),
                    "design":Counter({"Synchronic":2,"Allochronic":1}),
                    "experimental":Counter({"N":2,"Y":1}),
                    "record_count":3}}
        v=method_features(x,meta)
        self.assertEqual(v.shape,(1,10))
        self.assertAlmostEqual(v[0,0],2/3)
        self.assertAlmostEqual(v[0,1],1/3)
        self.assertAlmostEqual(v[0,2],np.log(3))
        self.assertAlmostEqual(v[0,3],np.log(4))
        self.assertEqual(v[0,9],0)

    def test_history_source_concentration_stable_and_nonnegative(self):
        rng=np.random.default_rng(20261114)
        X=rng.normal(size=(150,21));H=rng.normal(size=150)
        sp=np.array([f"S{k//5}" for k in range(150)])
        roles=np.array(["train" if k<75 else "heldout" for k in range(150)])
        q=source_history_concentration(H,X,sp,roles)
        self.assertGreater(q["history_unique_fraction_after_all_21_controls"],0)
        self.assertLessEqual(q["history_unique_fraction_after_all_21_controls"],1)
        self.assertGreater(q["effective_source_species_of_history_residual"],1)
        self.assertEqual(q["source_baseline_dim_including_intercept"],22)

    def test_original_known_truth_formula_on_synthetic_sites_only(self):
        rng=np.random.default_rng(20261114)
        X=rng.normal(size=(140,21));H=rng.normal(size=140)
        sp=np.array([f"S{k//5}" for k in range(140)])
        st=np.array([f"Study{k//5}" for k in range(140)])
        roles=np.array(["train" if k<70 else "heldout" for k in range(140)])
        contract=json.loads(CONTRACT.read_text())
        frozen=contract["synthetic_pre_response_qualification"]
        frozen["simulated_worlds_null"]=1
        frozen["simulated_worlds_positive"]=1
        worlds,passed=known_truth_gate(X,H,sp,st,roles,contract)
        self.assertEqual(worlds["null"]["worlds"],1)
        self.assertEqual(worlds["positive"]["worlds"],1)
        self.assertIsInstance(passed,bool)


if __name__=="__main__":
    unittest.main()
