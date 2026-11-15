"""Deterministic design tests; no biological outcomes or random simulation."""
import importlib.util
import json
import unittest
import numpy as np

SPEC = importlib.util.find_spec('history_design_bottleneck')
if SPEC is not None:
    import history_design_bottleneck as m

class ImplementationExists(unittest.TestCase):
    def test_core_exists(self):
        self.assertIsNotNone(SPEC, 'design bottleneck implementation missing')

@unittest.skipIf(SPEC is None, 'implementation is the red target')
class DesignTests(unittest.TestCase):
    def fixture(self):
        i=np.arange(96,dtype=float)
        X=np.column_stack([np.sin((k+1)*(i+.3)/17)+np.cos((k+2)*i/23) for k in range(21)])
        h=np.sin(i/3)+np.cos(i/7)+i/100
        sp=np.array(['s'+str(int(j)//12) for j in i])
        st=np.array(['p'+str(int(j)//6) for j in i])
        return X,h,sp,st
    def test_nested_energy_and_rank_budget(self):
        r=m.analyse_geometry(*self.fixture())
        s=r['spans']
        chain=['species_study','groups_environment11','groups_baseline21']
        for a,b in zip(chain,chain[1:]):
            self.assertLessEqual(s[b]['history_energy'],s[a]['history_energy']+1e-9)
            self.assertGreaterEqual(s[b]['rank'],s[a]['rank'])
        self.assertEqual(r['rank_budget']['residual_df_after_history'],96-s[chain[-1]]['rank']-1)
        self.assertFalse(r['biological_effect_tested'])
        json.dumps(r,allow_nan=False)
    def test_cell_projection_matches_direct_demeaning(self):
        X,h,sp,st=self.fixture()
        r=m.analyse_geometry(X,h,sp,st)
        direct=0
        for a,b in set(zip(sp,st)):
            z=h[(sp==a)&(st==b)]
            direct+=float(((z-z.mean())**2).sum())
        self.assertAlmostEqual(r['spans']['study_species_cells']['history_energy'],direct,places=9)
        self.assertEqual(r['cell_support']['occupied_cells'],16)
        self.assertEqual(r['cell_support']['repeated_cells'],16)
    def test_singleton_cells_have_no_within_cell_history(self):
        X,h,sp,st=self.fixture()
        st=np.array(['unique'+str(k) for k in range(96)])
        r=m.analyse_geometry(X,h,sp,st)
        self.assertEqual(r['cell_support']['singleton_cells'],96)
        self.assertEqual(r['cell_support']['cells_with_history_variation'],0)
        self.assertFalse(r['spans']['cells_baseline21']['history_identifiable'])
    def test_label_names_do_not_matter(self):
        X,h,sp,st=self.fixture()
        a=m.analyse_geometry(X,h,sp,st)
        b=m.analyse_geometry(X,h,['other'+s for s in sp],['other'+s for s in st])
        self.assertEqual(a['cell_support'],b['cell_support'])
        self.assertAlmostEqual(a['spans']['groups_baseline21']['history_energy'],b['spans']['groups_baseline21']['history_energy'],places=10)
    def test_duplicates_do_not_increase_span_information(self):
        X,h,_,_=self.fixture()
        a=m.span_summary(X,h)
        b=m.span_summary(np.column_stack([X,X[:,0],2*X[:,1]]),h)
        self.assertEqual(a['rank'],b['rank'])
        self.assertAlmostEqual(a['history_energy'],b['history_energy'],places=9)
    def test_bad_inputs_are_rejected(self):
        X,h,sp,st=self.fixture()
        for bad in [np.full_like(h,np.nan),h[:-1]]:
            with self.assertRaises(ValueError):m.analyse_geometry(X,bad,sp,st)
        with self.assertRaises(ValueError):m.analyse_geometry(X,h,sp[:-1],st)
        with self.assertRaises(ValueError):m.analyse_geometry(X,h,sp,['']*96)
    def test_no_information_from_a_group_constant(self):
        X,h,sp,st=self.fixture()
        h=np.repeat(np.arange(16,dtype=float),6)
        r=m.analyse_geometry(X,h,sp,st)
        self.assertFalse(r['spans']['species_study']['history_identifiable'])
    def test_matched_block_energy_formula(self):
        k=12;delta=1.7
        G=np.repeat(np.eye(k),2,axis=0)
        h=np.tile([-delta/2,delta/2],k)
        r=m.span_summary(G,h)
        self.assertAlmostEqual(r['history_energy'],k*delta**2/2,places=10)
        self.assertEqual(r['residual_df_after_history'],k-1)
    def test_minimum_blocks_is_actually_minimal(self):
        for d in (.5,1.,2.):
            r=m.required_blocks(d)
            self.assertGreaterEqual(r['power'],.8)
            self.assertLess(m.paired_power(r['blocks']-1,d),.8)
        self.assertGreater(m.required_blocks(.5)['blocks'],m.required_blocks(1)['blocks'])
    def test_power_null_equals_alpha_and_reject_bad_design(self):
        self.assertAlmostEqual(m.paired_power(20,1,gamma=0),1/60,places=12)
        for k,d in [(1,1),(20,0),(2.5,1)]:
            with self.assertRaises(ValueError):m.paired_power(k,d)
    def test_no_output_raw_labels_or_weights(self):
        X,h,sp,st=self.fixture()
        encoded=json.dumps(m.analyse_geometry(X,h,sp,st))
        self.assertNotIn('s0',encoded)
        self.assertNotIn('p0',encoded)
        self.assertNotIn('contrast_weights',encoded)

if __name__=='__main__': unittest.main()
