"""Synthetic-only first-shot Haldane numerical arithmetic & pre-response firewall tests."""
from __future__ import annotations
import math,sys,unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parent))
from proceed70_haldane_one_shot_empirical import (
    author_qualified_Haldane,only_after_authorized_response
)

def sample(**changes):
    row={"mean1":"10","mean2":"14","sd1":"2","sd2":"2",
         "n1":"10","n2":"10","generations":"2"}
    row.update(changes)
    return row


class HaldaneProspectiveOnlyTests(unittest.TestCase):
    def test_author_pooled_sd_and_per_generation_formula(self):
        a=author_qualified_Haldane(sample())
        self.assertAlmostEqual(a,math.log1p(1.),places=13)

    def test_guard_rejects_invalid_time_SD_and_study_sizes(self):
        for modify in (
            {"sd1":"0"},{"sd1":""},{"mean1":"0"},
            {"generations":"0"},{"generations":"nan"},{"n1":"1"},
            {"sd1":"-4"},{"mean2":"nan"},{"mean1":"1e9"},
            {"sd1":"1e-7"},{"mean2":"0"}
        ):
            self.assertIsNone(author_qualified_Haldane(sample(**modify)),modify)

    def test_exact_zero_change_rate_retained_not_cherrypicked(self):
        self.assertEqual(author_qualified_Haldane(sample(mean2="10")),0.)

    def test_frozen_response_opening_requires_upstream_pass(self):
        with self.assertRaisesRegex(ValueError,"No authorized"):
            only_after_authorized_response(b"mean1,mean2\n10,14\n",set(),
                                           {"status":"NOT_EVALUABLE_SYNTHETIC_RATE_PREDICTION_GATE_NO_RESPONSE_OPENING"})

    def test_negative_and_positive_mean_changes_have_same_abs_rate(self):
        a=author_qualified_Haldane(sample(mean1="10",mean2="14"))
        b=author_qualified_Haldane(sample(mean1="14",mean2="10"))
        self.assertAlmostEqual(a,b)


if __name__=="__main__":
    unittest.main()
