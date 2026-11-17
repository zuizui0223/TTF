"""Exact known-truth regression checks for working-model design variance."""
import copy,hashlib,json,sys,unittest
from fractions import Fraction
from pathlib import Path
HERE=Path(__file__).resolve().parent;sys.path.insert(0,str(HERE))
from ellis2007_quartz_design_precision_v01 import RECEIPT,SHA,evaluate,inverse,variance
class EllisPrecisionTests(unittest.TestCase):
    def setUp(self):self.data=json.loads(RECEIPT.read_text())
    def test_sha(self):self.assertEqual(hashlib.sha256(RECEIPT.read_bytes()).hexdigest(),SHA)
    def test_exact_inverse(self):
        self.assertEqual(inverse([[2,0],[0,3]]),[[Fraction(1,2),0],[0,Fraction(1,3)]])
    def test_singular_rejected(self):
        with self.assertRaisesRegex(ValueError,"singular"):inverse([[1,1],[1,1]])
    def test_one_arm_and_six_arm_variance(self):
        x=evaluate(self.data)
        self.assertAlmostEqual(x["by_restored_arm_count"]["1"]["min_variance_factor"],46/9)
        self.assertAlmostEqual(x["by_restored_arm_count"]["6"]["min_variance_factor"],2/3)
        self.assertAlmostEqual(x["one_arm_vs_all_six_variance_factor_ratio"],23/3)
    def test_exhaustive_subset_counts(self):
        x=evaluate(self.data)
        self.assertEqual([x["by_restored_arm_count"][str(k)]["subsets"] for k in range(1,7)],[6,15,20,15,6,1])
    def test_no_empirical_precision_claim(self):
        x=evaluate(self.data)
        self.assertFalse(x["original_fitness_genetics_GPS_LGM_outcomes_opened"])
        self.assertFalse(x["independent_replicate_or_causal_precision_claim"])
    def test_frozen_source_guard(self):
        d=copy.deepcopy(self.data);d["no_LGM_21ka_0BP_climate_effect_estimated"]=False
        with self.assertRaisesRegex(ValueError,"guard changed"):evaluate(d)
if __name__=="__main__":unittest.main()
