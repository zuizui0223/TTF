"""Seven source-only algebraic rank regression tests."""
import copy,json,sys,unittest
from pathlib import Path
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE))
from ellis2007_design_rank_v01 import rank,analyze,RECEIPT,SHA
import hashlib
class EllisDesignRankTests(unittest.TestCase):
    def setUp(self):
        self.data=json.loads(RECEIPT.read_text(encoding="utf-8"))
    def test_immutable_receipt(self):
        self.assertEqual(hashlib.sha256(RECEIPT.read_bytes()).hexdigest(),SHA)
    def test_exact_rank_and_scope(self):
        r=analyze(self.data)
        self.assertEqual((r["original_design_rank"],r["secondary_design_rank"]),(7,6))
        self.assertFalse(r["source_vs_garden_quartz_separately_identifiable_secondary"])
        self.assertFalse(r["fitness_genetics_climate_or_LGM_response_opened"])
    def test_known_truth_rational(self):
        self.assertEqual(rank([[1,0],[0,1],[1,1]]),2)
        self.assertEqual(rank([[1,1],[2,2]]),1)
    def test_duplicate_class_fails(self):
        d=copy.deepcopy(self.data)
        d["published_origin_site_pairs_by_class"][1]["class"]="local_same_source_and_garden"
        with self.assertRaisesRegex(ValueError,"duplicate"):analyze(d)
    def test_class_drift_fails(self):
        d=copy.deepcopy(self.data)
        d["published_origin_site_pairs_by_class"][1]["original_Hargreaves_extract_observed"]=1
        with self.assertRaisesRegex(ValueError,"topology drift"):analyze(d)
    def test_no_outcome_guards(self):
        for key in ("no_fitness_growth_survival_reproduction_or_AFLP_response_rows_read",
                    "no_LGM_21ka_0BP_climate_effect_estimated",
                    "old_Bontrager_44of50_terminal_stop_unchanged"):
            d=copy.deepcopy(self.data);d[key]=False
            with self.assertRaisesRegex(ValueError,"guard changed"):analyze(d)
    def test_ragged_matrix_fails(self):
        with self.assertRaisesRegex(ValueError,"ragged"):rank([[1,2],[1]])
if __name__=="__main__":unittest.main()
