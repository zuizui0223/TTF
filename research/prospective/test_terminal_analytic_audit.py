import unittest
import numpy as np
import terminal_analytic_audit as a

class AnalyticAuditTests(unittest.TestCase):
    def test_wilson_frozen_failure(self):
        lo, hi = a.wilson(7,100)
        self.assertAlmostEqual(lo,0.03431926106727268)
        self.assertAlmostEqual(hi,0.137495147390735)
        self.assertEqual(max(k for k in range(101) if a.wilson(k,100)[1] <= .1),4)

    def test_duplicate_feature_has_no_new_information_but_changes_ridge(self):
        x=np.array([-1.,1.,-1.,1.])[:,None]
        train=np.array([True,True,False,False]);w=np.array([.5,.5])
        p=a.ridge_operator(x,train,w)
        q=a.ridge_operator(np.column_stack([x,x]),train,w)
        target=np.array([-1.,1.]); original=target
        e0=np.mean((target-p@original)**2);e1=np.mean((target-q@original)**2)
        self.assertAlmostEqual(1-e1/e0,0.7256235827664399)

    def test_operator_agrees_with_independent_direct_fit(self):
        X=np.column_stack([np.arange(12),np.sin(np.arange(12))])
        train=np.arange(12)<8; w=np.arange(1,9,dtype=float);w/=w.sum()
        y=np.cos(np.arange(8))
        p=a.ridge_operator(X,train,w)
        mu=w@X[train];sd=np.sqrt(w@((X[train]-mu)**2));z=(X-mu)/sd
        ym=w@y
        b=np.linalg.solve(z[train].T@(z[train]*w[:,None])+.1*np.eye(2),z[train].T@(w*(y-ym)))
        np.testing.assert_allclose(p@y,ym+z[~train]@b,atol=1e-12)

    def test_gaussian_square_covariance_scalar(self):
        val=a.square_covariance(np.array([2.]),np.array([3.]),np.array([[.5]]))
        self.assertAlmostEqual(val[0,0],12.5)

    def test_shared_study_cross_species_covariance_is_retained(self):
        v=a.noise_covariance(np.array(['a','b','c']),np.array(['s','s','t']))
        np.testing.assert_allclose(np.diag(v),.75)
        self.assertEqual(v[0,1],.25)
        self.assertEqual(v[0,2],0)

    def test_identical_operators_have_zero_paired_loss_covariance(self):
        A=np.eye(3);m=np.array([0.,1.,2.]);v=np.eye(3)*.75
        mean,cov=a.loss_difference_moments(A,A,m,v)
        np.testing.assert_allclose(mean,0,atol=1e-12)
        np.testing.assert_allclose(cov,0,atol=1e-12)

    def test_error_operator_matches_residual(self):
        train=np.array([True,False,True,False])
        p=np.array([[.3,.7],[.6,.4]])
        y=np.array([4.,3.,1.,2.])
        A=a.error_operator(p,train)
        np.testing.assert_allclose(A@y,y[~train]-p@y[train])

    def test_no_implementation_random_draws(self):
        import inspect
        self.assertNotIn('np.random',inspect.getsource(a))

if __name__=='__main__': unittest.main()
