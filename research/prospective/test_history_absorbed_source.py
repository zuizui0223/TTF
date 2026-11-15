"""Source guards and aggregate-only reporting, using deterministic invented X/H."""
import ast
import inspect
import unittest
import numpy as np
try:
    import history_absorbed_source as runner
except ModuleNotFoundError:
    runner=None


def source_fixture():
    sp=[];st=[];roles=[]
    for arm,sizes,split_at in [('train',[4]*4+[3]*25,19),('heldout',[3]*20+[2]*7,18)]:
        for j,m in enumerate(sizes):
            for k in range(m):
                sp.append(f'{arm}-sp{j}')
                st.append(f'{arm}-st{j}-'+('1' if j<split_at and k>0 else '0'))
                roles.append(arm)
    i=np.arange(165)[:,None]+1
    j=np.arange(21)[None,:]+1
    X=np.sin(i*j*.197)+.3*np.cos(i*(j+1)*.317)
    H=np.sin(np.arange(165)*1.721)+np.arange(165)*.01
    return X,H,np.asarray(sp),np.asarray(st),np.asarray(roles)

class RunnerExists(unittest.TestCase):
    def test_source_runner_available(self):
        self.assertIsNotNone(runner,'unknown-scale source runner missing')

@unittest.skipIf(runner is None,'source runner not implemented')
class SourceTests(unittest.TestCase):
    def test_runtime_selection_requires_both_original_hashes(self):
        self.assertTrue(hasattr(runner,'select_exact_mode'),'exact-mode selector is missing')
        reports={'wrong':{'X_sha256':'a','H_sha256':'H'},
                 'exact':{'X_sha256':'X','H_sha256':'H'}}
        self.assertEqual(runner.select_exact_mode(reports,'X','H'),'exact')
        self.assertIsNone(runner.select_exact_mode(reports,'Z','H'))
        self.assertIsNone(runner.select_exact_mode({'failed':{'worker_failed':True}},'X','H'))

    def test_source_aggregate_has_no_raw_arrays(self):
        r=runner.analyse_source(*source_fixture())
        self.assertEqual(r['heldout_rows'],74)
        self.assertEqual(r['heldout_species'],27)
        self.assertEqual(r['heldout_studies'],45)
        self.assertFalse(r['real_phenotype_validity_established'])
        self.assertFalse(r['empirical_opening_authorized'])
        def noarrays(obj):
            self.assertFalse(isinstance(obj,np.ndarray))
            if isinstance(obj,dict):
                for v in obj.values():noarrays(v)
            elif isinstance(obj,list):
                for v in obj:noarrays(v)
        noarrays(r)
        self.assertNotIn('unit_contrast',r['within_group_design'])
        self.assertNotIn('residual_projector',r['within_group_design'])

    def test_shape_or_finiteness_change_stops(self):
        X,H,s,st,role=source_fixture()
        with self.assertRaises(ValueError):runner.analyse_source(X[:-1],H,s,st,role)
        H[4]=np.nan
        with self.assertRaises(ValueError):runner.analyse_source(X,H,s,st,role)

    def test_exact_original_source_hash_guards(self):
        a=np.array([1.,2.])
        self.assertTrue(runner.matches_original(a,a,runner.digest(a),runner.digest(a)))
        b=a.copy();b[1]=np.nextafter(b[1],np.inf)
        self.assertFalse(runner.matches_original(a,b,runner.digest(a),runner.digest(a)))

    def test_no_outcome_or_simulator_imports(self):
        s=inspect.getsource(runner)
        tree=ast.parse(s)
        funcs={node.func.id for node in ast.walk(tree) if isinstance(node,ast.Call) and isinstance(node.func,ast.Name)}
        self.assertFalse(funcs & {'known_truth_labels','known_truth_gate','ridge_predict','run_first_shot'})
        self.assertNotIn('np.random',s)

if __name__=='__main__':unittest.main()
