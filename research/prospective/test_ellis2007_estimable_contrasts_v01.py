"""Source-only Ellis estimable-contrast known-truth checks."""
import copy,hashlib,json,sys,unittest
from pathlib import Path
HERE=Path(__file__).resolve().parent;sys.path.insert(0,str(HERE))
from ellis2007_estimable_contrasts_v01 import RECEIPT,SHA,evaluate,null_basis,estimable
class ContrastTests(unittest.TestCase):
    def setUp(self):self.data=json.loads(RECEIPT.read_text())
    def test_sha(self):self.assertEqual(hashlib.sha256(RECEIPT.read_bytes()).hexdigest(),SHA)
    def test_rank_and_null_vector(self):
        a=evaluate(self.data)
        self.assertEqual((a["primary_rank"],a["secondary_rank"]),(7,6))
        self.assertEqual(a["secondary_null_vector"],[0,-1,1,0,0,0,0])
    def test_individual_quartz_unidentified(self):
        a=evaluate(self.data)["contrasts"]
        for k in ("source_quartz_only","garden_quartz_only","source_minus_garden_quartz"):
            self.assertTrue(a[k]["primary_estimable"])
            self.assertFalse(a[k]["secondary_estimable"])
    def test_sum_and_basin_remain_estimable(self):
        a=evaluate(self.data)["contrasts"]
        for k in ("source_plus_garden_quartz","source_basin_2_vs_1","garden_basin_2_vs_1"):
            self.assertTrue(a[k]["secondary_estimable"])
    def test_known_truth_null_basis(self):
        self.assertEqual(null_basis([[1,1],[2,2]]),[[-1,1]])
        self.assertTrue(estimable([1,1],[[-1,1]]))
        self.assertFalse(estimable([1,0],[[-1,1]]))
    def test_source_only_protections(self):
        a=evaluate(self.data)
        self.assertTrue(a["all_original_outcomes_unopened"])
        self.assertTrue(a["frozen_endpoints_unchanged"])
    def test_invalid_source_gate_fails(self):
        d=copy.deepcopy(self.data)
        d["no_LGM_21ka_0BP_climate_effect_estimated"]=False
        with self.assertRaisesRegex(ValueError,"guard changed"):evaluate(d)
if __name__=="__main__":unittest.main()
