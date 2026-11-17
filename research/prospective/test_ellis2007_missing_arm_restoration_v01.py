"""Seven tests for six-arm restoration enumeration, with fail-closed source-only gates."""
import copy, hashlib, json, sys, unittest
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from ellis2007_missing_arm_restoration_v01 import RECEIPT, SHA, evaluate

class MissingArmRestorationTests(unittest.TestCase):
    def setUp(self):
        self.data = json.loads(RECEIPT.read_text(encoding="utf-8"))

    def test_immutable_input(self):
        self.assertEqual(hashlib.sha256(RECEIPT.read_bytes()).hexdigest(), SHA)

    def test_all_sixty_four_subsets(self):
        result = evaluate(self.data)
        self.assertEqual(result["subset_count"], 64)
        self.assertEqual([result["by_number_of_restored_arms"][str(k)]["subsets"] for k in range(7)],
                         [1, 6, 15, 20, 15, 6, 1])

    def test_one_arm_restores_algebraic_rank_and_connectivity(self):
        result = evaluate(self.data)
        self.assertEqual(result["minimum_added_original_arms_for_full_rank"], 1)
        self.assertEqual(result["minimum_added_original_arms_for_connected_graph"], 1)
        self.assertEqual(result["by_number_of_restored_arms"]["1"]["full_rank_subsets"], 6)
        self.assertEqual(result["by_number_of_restored_arms"]["1"]["connected_subsets"], 6)

    def test_no_restoration_does_not_repair_rank(self):
        result = evaluate(self.data)
        self.assertEqual(result["by_number_of_restored_arms"]["0"]["full_rank_subsets"], 0)
        self.assertEqual(result["baseline_components"], 2)

    def test_no_biological_claim(self):
        result = evaluate(self.data)
        self.assertFalse(result["actual_primary_response_values_recovered"])
        self.assertFalse(result["restored_one_arm_sufficient_for_precision_or_causal_identification"])
        self.assertFalse(result["fitness_genetics_LGM_or_adaptation_response_opened"])
        self.assertTrue(result["frozen_TTF_endpoints_unchanged"])

    def test_class_count_drift_fails_closed(self):
        data = copy.deepcopy(self.data)
        data["published_origin_site_pairs_by_class"][1]["original_Hargreaves_extract_observed"] = 1
        with self.assertRaisesRegex(ValueError, "topology drift"):
            evaluate(data)

    def test_outcome_guard_fails_closed(self):
        data = copy.deepcopy(self.data)
        data["no_fitness_growth_survival_reproduction_or_AFLP_response_rows_read"] = False
        with self.assertRaisesRegex(ValueError, "guard changed"):
            evaluate(data)

if __name__ == "__main__":
    unittest.main()
