"""Source boundary tests; deterministic artificial inputs only."""
import ast
import hashlib
import importlib.util
import json
from pathlib import Path
import unittest
import numpy as np

SPEC=importlib.util.find_spec('history_design_bottleneck_source')
if SPEC is not None:
    import history_design_bottleneck_source as m

class RunnerExists(unittest.TestCase):
    def test_runner_exists(self):
        self.assertIsNotNone(SPEC,'new source-anatomy runner missing')

@unittest.skipIf(SPEC is None,'source runner missing is the red target')
class BoundaryTests(unittest.TestCase):
    def fixture(self):
        return np.arange(165*21,dtype=float).reshape(165,21),np.arange(165,dtype=float)
    def test_exact_original_hashes_both_required(self):
        X,H=self.fixture();a=m.digest(X);b=m.digest(H)
        self.assertTrue(m.require_identity(X,H,a,b))
        for name in ('X','H'):
            x=X.copy();h=H.copy()
            target=x if name=='X' else h
            target.flat[0]=np.nextafter(target.flat[0],1.)
            with self.assertRaises(ValueError):m.require_identity(x,h,a,b)
    def test_wrong_shape_or_nonfinite_abstains(self):
        X,H=self.fixture()
        with self.assertRaises(ValueError):m.require_identity(X[:-1],H,m.digest(X),m.digest(H))
        X[0,0]=np.nan
        with self.assertRaises(ValueError):m.require_identity(X,H,m.digest(X),m.digest(H))
    def test_no_simulation_or_real_response_function_calls(self):
        root=ast.parse(Path(m.__file__).read_text())
        forbidden={'known_truth_gate','known_truth_labels','ridge_predict','missingness_summary','random','rvs'}
        calls=set()
        for x in ast.walk(root):
            if isinstance(x,ast.Call):
                calls.add(getattr(x.func,'id',getattr(x.func,'attr','')))
        self.assertFalse(calls & forbidden)
    def test_terminal_and_interpretation_guards_are_present(self):
        text=Path(m.__file__).read_text()
        self.assertIn('before!=terminal_path.read_bytes()',text)
        self.assertIn("'empirical_opening_authorized':False",text)
        self.assertIn("'new_random_worlds':0",text)

if __name__=='__main__':unittest.main()
