"""Stage23 synthetic-only mathematical tests. No Jensen source/response loading."""
import unittest
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import numpy as np
from stage23_paired_case_climate_shift import case_history_magnitude,pair_summary

class CaseClimateTests(unittest.TestCase):
    def test_locked_four_band_scaled_magnitude(self):
        x=np.array([[1,2,3,4,0,0,0,0],[0,0,0,0,0,0,0,0]],float)
        np.testing.assert_allclose(case_history_magnitude(x,[1,2,3,4]),[2.,0.])
        with self.assertRaises(ValueError):
            case_history_magnitude(x,[1,-1,3,4])
        with self.assertRaises(ValueError):
            case_history_magnitude(x[:,:7],[1,2,3,4])

    def test_identical_jitter_means_zero_predictor_change(self):
        x=np.zeros((12,8),float)
        x[:,0]=np.linspace(0,11,12)
        x[:,4:8]=1
        records=[{"scientific_name":"Testus syntheticus","case_points":12,"split_role":0}]
        out=pair_summary(x,x,np.zeros(12,dtype=int),records,
                         ["bio01_LGM","bio07_LGM","bio12_LGM","bio15_LGM",
                          "bio01_current","bio07_current","bio12_current","bio15_current"],
                         [1,1,1,1])
        self.assertEqual(out["original_case_records"],12)
        self.assertEqual(out["paired_abs_history_shift_median"],0.)
        self.assertEqual(out["paired_abs_history_shift_p95"],0.)
        self.assertEqual(out["original_moved_history_spearman"],1.)
        self.assertEqual(out["species_with_undefined_spearman"],0)
        self.assertEqual(out["per_species_aggregates"][0]["paired_abs_history_shift_median"],0.)

    def test_different_case_site_count_cannot_be_used(self):
        a=np.zeros((4,8))
        with self.assertRaisesRegex(ValueError,"Frozen case count"):
            pair_summary(a,a,np.zeros(4,dtype=int),
                         [{"scientific_name":"Testus syntheticus","case_points":5,"split_role":0}],
                         [str(j) for j in range(8)],[1,1,1,1])

    def test_missingness_abstains(self):
        a=np.zeros((12,8))
        a[1,2]=np.nan
        with self.assertRaisesRegex(ValueError,"No imputation"):
            pair_summary(a,a,np.zeros(12,dtype=int),
                         [{"scientific_name":"Testus syntheticus","case_points":12,"split_role":0}],
                         [str(j) for j in range(8)],[1,1,1,1])

if __name__=="__main__":
    unittest.main()
