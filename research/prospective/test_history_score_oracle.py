"""Deterministic mathematical controls only; zero Monte Carlo or real outcomes."""
import ast
import importlib.util
import sys
import unittest
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
AVAILABLE = importlib.util.find_spec('history_score_oracle') is not None
if AVAILABLE:
    import history_score_oracle as core

class ImplementationExists(unittest.TestCase):
    def test_oracle_implementation_available(self):
        self.assertTrue(AVAILABLE, 'history_score_oracle source-only contrast not yet implemented')

@unittest.skipUnless(AVAILABLE, 'implementation not yet present')
class OracleTests(unittest.TestCase):
    def setUp(self):
        t=np.arange(24.,dtype=float)
        self.Z=np.column_stack([np.ones(24),t,np.sin(t/3)])
        self.h=np.cos(t/2)+t/20
        self.species=np.repeat(np.arange(8),3)
        self.studies=np.repeat(np.arange(4),6)
        self.S=core.cluster_covariance(self.species,self.studies)

    def test_nuisance_annihilation_and_unit_variance(self):
        obj=core.history_contrast(self.Z,self.h,self.S)
        np.testing.assert_allclose(obj['z_weights']@self.Z,0,atol=1e-11)
        u=obj['z_weights']
        self.assertAlmostEqual(float(u@self.S@u),1,places=12)
        self.assertAlmostEqual(float(u@self.h),np.sqrt(obj['information']),places=11)

    def test_gls_schur_complement_is_independent_check(self):
        obj=core.history_contrast(self.Z,self.h,self.S)
        P=np.linalg.inv(self.S)
        R=P-P@self.Z@np.linalg.solve(self.Z.T@P@self.Z,self.Z.T@P)
        direct=float(self.h@R@self.h)
        self.assertAlmostEqual(obj['information'],direct,places=10)

    def test_duplicate_nuisance_does_not_change_information(self):
        a=core.history_contrast(self.Z,self.h,self.S)
        dup=np.column_stack([self.Z,self.Z[:,1],2*self.Z[:,2]])
        b=core.history_contrast(dup,self.h,self.S)
        np.testing.assert_allclose(a['z_weights'],b['z_weights'],rtol=1e-10,atol=1e-11)
        self.assertEqual(a['nuisance_rank'],b['nuisance_rank'])

    def test_duplicate_history_is_nonidentifiable_not_false_gain(self):
        with self.assertRaisesRegex(ValueError,'not identifiable'):
            core.history_contrast(self.Z,3*self.Z[:,1]-self.Z[:,2],self.S)

    def test_reparameterizing_nuisance_span_changes_nothing(self):
        a=core.history_contrast(self.Z,self.h,self.S)
        T=np.array([[2,1,0],[0,0.1,0.4],[0,0,3.]])
        b=core.history_contrast(self.Z@T,self.h,self.S)
        np.testing.assert_allclose(a['z_weights'],b['z_weights'],rtol=1e-10,atol=1e-11)

    def test_covariance_scaling_has_known_information_effect(self):
        a=core.history_contrast(self.Z,self.h,self.S)
        b=core.history_contrast(self.Z,self.h,4*self.S)
        self.assertAlmostEqual(a['information']/4,b['information'],places=11)

    def test_history_scaling_and_nuisance_addition(self):
        a=core.history_contrast(self.Z,self.h,self.S)
        b=core.history_contrast(self.Z,2*self.h+self.Z@np.array([5,-1,2]),self.S)
        self.assertAlmostEqual(4*a['information'],b['information'],places=9)
        np.testing.assert_allclose(a['z_weights'],b['z_weights'],atol=1e-11)

    def test_cluster_covariance_includes_shared_study(self):
        self.assertEqual(self.S[0,1],0.5)
        self.assertEqual(self.S[0,4],0.25)
        self.assertEqual(self.S[0,8],0)
        self.assertEqual(self.S[0,0],0.75)

    def test_exact_null_size_and_monotonic_power(self):
        alpha=1/60
        self.assertAlmostEqual(core.two_sided_power(0,1,alpha),alpha,places=14)
        vals=[core.two_sided_power(d,1,alpha) for d in [0,.5,1,2,3]]
        self.assertTrue(np.all(np.diff(vals)>0))
        self.assertAlmostEqual(core.two_sided_power(-2,1,alpha),vals[3],places=14)
        d=core.effect_for_power(4.0,0.8,alpha)
        self.assertAlmostEqual(core.two_sided_power(d*2,1,alpha),.8,places=11)

    def test_wrong_covariance_size_is_not_forced_to_alpha(self):
        a=core.history_contrast(self.Z,self.h,np.eye(len(self.h)))
        actual=float(a['z_weights']@self.S@a['z_weights'])
        rate=core.two_sided_power(0,actual,1/60)
        self.assertGreater(abs(rate-1/60),1e-5)

    def test_bounded_shape_and_covariance_validation(self):
        for cov in [-self.S,self.S[:10,:10],np.full_like(self.S,np.nan)]:
            with self.assertRaises(ValueError): core.history_contrast(self.Z,self.h,cov)
        asym=self.S.copy();asym[0,1]+=0.1
        with self.assertRaises(ValueError):core.history_contrast(self.Z,self.h,asym)
        with self.assertRaises(ValueError):core.history_contrast(self.Z,self.h[:10],self.S)
        bad=self.h.copy();bad[0]=np.nan
        with self.assertRaises(ValueError):core.history_contrast(self.Z,bad,self.S)
        with self.assertRaises(ValueError):core.two_sided_power(0,-1,.05)
        with self.assertRaises(ValueError):core.effect_for_power(0,.8,.05)

    def test_core_has_no_response_or_rng_interface(self):
        import inspect
        self.assertEqual(list(inspect.signature(core.history_contrast).parameters),['nuisance','history','covariance'])
        tree=ast.parse(Path(core.__file__).read_text())
        names=[]
        for node in ast.walk(tree):
            if isinstance(node,ast.Attribute): names.append(node.attr)
        self.assertNotIn('default_rng',names)
        self.assertNotIn('rvs',names)
        self.assertNotIn('random',names)

if __name__=='__main__':unittest.main()
