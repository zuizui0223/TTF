"""Deterministic source-design validation, no original response or simulations."""
import ast
import importlib.util
from pathlib import Path
import sys
import unittest
import numpy as np
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE))
AVAILABLE=importlib.util.find_spec('history_score_oracle_source') is not None
if AVAILABLE: import history_score_oracle_source as run

class RunnerExists(unittest.TestCase):
    def test_source_design_runner_exists(self):
        self.assertTrue(AVAILABLE,'source design oracle runner not implemented')

@unittest.skipUnless(AVAILABLE,'source design runner not implemented')
class SourceTests(unittest.TestCase):
    def fixture(self):
        # Fixed trigonometric source features, not simulated outcomes.
        t=np.arange(165.)
        X=np.column_stack([np.sin((t+1)*(j+1)/23)+np.cos((t+2)*(j+2)/41) for j in range(21)])
        H=np.sin(t/2)+np.cos(t/5)+4
        role=np.array(['train']*91+['heldout']*74)
        sp=np.concatenate([np.array([f'a{i%29}' for i in range(91)]),np.array([f'b{i%27}' for i in range(74)])])
        st=np.concatenate([np.array([f'sa{i%48}' for i in range(91)]),np.array([f'sb{i%45}' for i in range(74)])])
        return X,H,sp,st,role

    def test_report_is_aggregate_only_and_unconditional_authorization_false(self):
        d=run.analyse_design(*self.fixture())
        self.assertEqual(d['heldout_rows'],74)
        self.assertEqual(d['heldout_species'],27)
        self.assertFalse(d['empirical_response_opening_authorized'])
        self.assertFalse(d['oracle_is_real_biological_validation'])
        self.assertNotIn('z_weights',str(d))
        self.assertAlmostEqual(d['known_covariance_oracle']['analytic_null_rejection_probability'],1/60,places=12)
        self.assertEqual(d['known_covariance_oracle']['power_grid'][0]['gamma_per_train_history_SD'],0)

    def test_original_source_shape_drift_stops(self):
        X,H,s,st,r=self.fixture()
        with self.assertRaises(ValueError):run.analyse_design(X[:,:20],H,s,st,r)
        bad=r.copy();bad[0]='heldout'
        with self.assertRaises(ValueError):run.analyse_design(X,H,s,st,bad)
        st[0]=st[-1]
        with self.assertRaisesRegex(ValueError,'leakage'):run.analyse_design(X,H,s,st,r)

    def test_exact_source_hash_no_tolerance_for_one_ulp(self):
        a=np.array([1.,2.]);b=a.copy();b[0]=np.nextafter(b[0],np.inf)
        self.assertNotEqual(run.digest(a),run.digest(b))

    def test_runner_does_not_call_old_simulations_or_outcome_opening(self):
        tree=ast.parse(Path(run.__file__).read_text())
        calls=[]
        for node in ast.walk(tree):
            if isinstance(node,ast.Call):
                fn=node.func
                calls.append(fn.id if isinstance(fn,ast.Name) else fn.attr if isinstance(fn,ast.Attribute) else '')
        self.assertFalse(set(calls)&{'known_truth_gate','known_truth_labels','ridge_predict','default_rng','rvs'})

if __name__=='__main__':unittest.main()
