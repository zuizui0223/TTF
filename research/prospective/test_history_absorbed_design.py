"""Deterministic design mathematics only; no source responses or random draws."""
import inspect
import unittest
import numpy as np
try:
    import history_absorbed_design as core
except ModuleNotFoundError:
    core = None


def fixture():
    x = np.arange(80, dtype=float)
    species = np.array(['sp'+str(i//8) for i in range(80)])
    studies = np.array(['st'+str(i//2) for i in range(80)])
    Z = np.column_stack([np.ones(80), np.sin(x*.37), np.cos(x*.81)])
    h = np.sin(x*1.217)+.7*np.cos(x*.93)
    return Z,h,species,studies


class ImplementationExists(unittest.TestCase):
    def test_design_implementation_available(self):
        self.assertIsNotNone(core, 'unknown-variance absorption design implementation is missing')


@unittest.skipIf(core is None, 'implementation not yet written')
class AbsorbedDesignTests(unittest.TestCase):
    def setUp(self):
        self.Z,self.h,self.s,self.st = fixture()
        self.d = core.absorbed_design(self.Z,self.h,self.s,self.st)

    def test_rank_and_residual_degrees(self):
        self.assertEqual(self.d['status'],'IDENTIFIABLE_UNKNOWN_SCALE_GAUSSIAN_DESIGN')
        self.assertEqual(self.d['group_incidence_rank'],40)
        self.assertEqual(self.d['nuisance_rank'],42)
        self.assertEqual(self.d['residual_df'],37)
        self.assertGreater(self.d['history_residual_energy'],0)

    def test_annihilation_and_projection(self):
        d=self.d; u=d['unit_contrast']; M=d['residual_projector']; W=d['nuisance_matrix']
        np.testing.assert_allclose(u@W,0,atol=2e-12)
        np.testing.assert_allclose(M@W,0,atol=3e-12)
        np.testing.assert_allclose(M@self.h,0,atol=3e-12)
        np.testing.assert_allclose(M@M,M,atol=3e-12)
        self.assertAlmostEqual(float(u@u),1.,12)
        self.assertAlmostEqual(float(np.trace(M)),37.,10)

    def test_arbitrary_unknown_group_variances_cancel(self):
        A=core.incidence(self.s); B=core.incidence(self.st)
        a=np.arange(A.shape[1],dtype=float)+1
        b=(np.arange(B.shape[1],dtype=float)+1)**2
        S=.09*np.eye(80)+(A*a)@A.T+(B*b)@B.T
        u=self.d['unit_contrast'];M=self.d['residual_projector']
        self.assertAlmostEqual(float(u@S@u),.09,9)
        np.testing.assert_allclose(M@S@M,.09*M,atol=2e-10)
        np.testing.assert_allclose(u@S@M,0,atol=2e-10)

    def test_fixed_offsets_and_unknown_scale_cancel_in_pivot(self):
        d=self.d;u=d['unit_contrast'];M=d['residual_projector'];df=d['residual_df']
        y=np.sin(np.arange(80)*1.71)+.2*self.h
        shift=core.incidence(self.s)@np.arange(10)+core.incidence(self.st)@np.arange(40)
        pivot=lambda v:float(u@v)/np.sqrt(float(v@M@v)/df)
        self.assertAlmostEqual(pivot(y),pivot(7*y+shift),9)

    def test_duplicates_and_reparametrization_no_new_information(self):
        dup=np.column_stack([self.Z,self.Z[:,1],self.Z[:,2]*7])
        d=core.absorbed_design(dup,self.h,self.s,self.st)
        self.assertAlmostEqual(d['history_residual_energy'],self.d['history_residual_energy'],10)
        self.assertEqual(d['nuisance_rank'],self.d['nuisance_rank'])
        np.testing.assert_allclose(d['residual_projector'],self.d['residual_projector'],atol=2e-12)

    def test_history_in_group_span_is_not_identifiable(self):
        h=core.incidence(self.s)@np.arange(10)+self.Z[:,1]
        d=core.absorbed_design(self.Z,h,self.s,self.st)
        self.assertEqual(d['status'],'NO_WITHIN_GROUP_HISTORY_INFORMATION')
        self.assertNotIn('unit_contrast',d)

    def test_saturated_model_does_not_emit_a_t_test(self):
        x=np.arange(8,dtype=float)
        d=core.absorbed_design(np.ones((8,1)),x,[str(i) for i in range(8)],['all']*8)
        self.assertNotEqual(d['status'],'IDENTIFIABLE_UNKNOWN_SCALE_GAUSSIAN_DESIGN')
        self.assertNotIn('unit_contrast',d)

    def test_alpha_and_power_are_finite_and_monotone(self):
        energy=self.d['history_residual_energy'];df=self.d['residual_df']
        powers=[core.student_power(g,energy,df,.5) for g in [0,.25,.5,1]]
        self.assertAlmostEqual(powers[0],1/60,12)
        self.assertTrue(all(a<b for a,b in zip(powers,powers[1:])))
        g=core.detectable_effect(energy,df,.5)
        self.assertAlmostEqual(core.student_power(g,energy,df,.5),.8,9)
        self.assertAlmostEqual(core.detectable_effect(energy,df,1.),2*g,10)

    def test_bad_shapes_values_and_labels_fail(self):
        for Z,h,s,st in [(self.Z[:-1],self.h,self.s,self.st),
                         (self.Z,self.h[:-1],self.s,self.st),
                         (self.Z,np.full(80,np.nan),self.s,self.st),
                         (self.Z,self.h,self.s,['']*80)]:
            with self.assertRaises(ValueError):core.absorbed_design(Z,h,s,st)
        for pars in [(0,10,.5),(1,0,.5),(1,5,0),(float('inf'),5,.5)]:
            with self.assertRaises(ValueError):core.student_power(.5,*pars)

    def test_group_label_names_do_not_change_design(self):
        d=core.absorbed_design(self.Z,self.h,['renamed_'+s for s in self.s],['z'+s for s in self.st])
        self.assertAlmostEqual(d['history_residual_energy'],self.d['history_residual_energy'],10)

    def test_no_rng_or_biological_response_interface(self):
        src=inspect.getsource(core)
        self.assertNotIn('np.random',src)
        self.assertNotIn('read_csv',src)
        sig=inspect.signature(core.absorbed_design)
        self.assertEqual(list(sig.parameters),['Z','h','species','studies'])

if __name__=='__main__':unittest.main()
