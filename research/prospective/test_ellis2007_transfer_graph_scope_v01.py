"""Regression checks for a frozen source-only Ellis transfer graph diagnosis."""
from __future__ import annotations
import copy
import hashlib
import json
import sys
import unittest
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from ellis2007_transfer_graph_scope_v01 import RECEIPT, SHA256, analyze


class Ellis2007TransferGraphTests(unittest.TestCase):
    def setUp(self):
        self.data = json.loads(RECEIPT.read_text(encoding="utf-8"))

    def test_original_immutable_receipt(self):
        self.assertEqual(hashlib.sha256(RECEIPT.read_bytes()).hexdigest(), SHA256)

    def test_comparison_graph_exact_known_truth(self):
        r = analyze(self.data)
        self.assertEqual((r["original_directed_arms"], r["secondary_directed_arms"]), (24, 18))
        self.assertEqual((r["original_connected_components"], r["secondary_connected_components"]), (1, 2))
        self.assertEqual(r["secondary_component_sizes"], [3, 3])
        self.assertTrue(r["secondary_components_homogeneous_by_quartz"])
        self.assertEqual((r["original_within_basin_cross_quartz_arms"],
                          r["secondary_within_basin_cross_quartz_arms"]), (6, 0))
        self.assertEqual((r["original_observed_arm_classes"], r["secondary_observed_arm_classes"]), (3, 2))
        self.assertFalse(r["secondary_within_basin_cross_quartz_comparison_present"])
        self.assertFalse(r["genetic_or_fitness_or_LGM_response_opened"])
        self.assertFalse(r["former_frozen_endpoints_changed"])
        self.assertFalse(r["source_garden_IDs_or_coordinates_exported"])

    def test_partial_category_fails_closed_without_guessing_ids(self):
        data = copy.deepcopy(self.data)
        row = data["published_origin_site_pairs_by_class"][1]
        row["original_Hargreaves_extract_observed"] = 1
        row["not_represented_in_secondary_original"] = 5
        data["observed_Hargreaves_distinct_directed_source_garden_pairs"] = 19
        with self.assertRaisesRegex(ValueError, "partial class"):
            analyze(data)

    def test_duplicate_class_rejected(self):
        data = copy.deepcopy(self.data)
        data["published_origin_site_pairs_by_class"][1]["class"] = "local_same_source_and_garden"
        with self.assertRaisesRegex(ValueError, "duplicate or unknown"):
            analyze(data)

    def test_primary_expected_counts_immutable(self):
        data = copy.deepcopy(self.data)
        row = data["published_origin_site_pairs_by_class"][1]
        row["paper_expected"] = 5
        row["not_represented_in_secondary_original"] = 5
        data["original_primary_expected_directed_source_garden_pairs"] = 23
        with self.assertRaisesRegex(ValueError, "primary design drift"):
            analyze(data)

    def test_no_response_read_guard_required(self):
        data = copy.deepcopy(self.data)
        data["no_fitness_growth_survival_reproduction_or_AFLP_response_rows_read"] = False
        with self.assertRaisesRegex(ValueError, "frozen protection changed"):
            analyze(data)

    def test_no_lgm_or_previous_stop_guard_required(self):
        for key in ("no_LGM_21ka_0BP_climate_effect_estimated",
                    "old_Bontrager_44of50_terminal_stop_unchanged"):
            data = copy.deepcopy(self.data)
            data[key] = False
            with self.assertRaisesRegex(ValueError, "frozen protection changed"):
                analyze(data)


if __name__ == "__main__":
    unittest.main()
