"""Post-outcome mathematical diagnostics, not empirical hypothesis tests.

Only synthetic numerical arrays and the public frozen terminal receipt are
read. These checks do not replace the original 259-species estimand, open
any additional nucleotide data, or summarize the private 255 coefficients.
"""
from pathlib import Path
import hashlib
import json
import unittest

import numpy as np


class ResponseScaleMathematics(unittest.TestCase):
    def setUp(self):
        self.h = np.array([-1., -1., -1., -1., 1., 1., 1., 1.])
        self.s = np.array([-1., -1., 1., 1., -1., -1., 1., 1.])
        self.c = np.array([-1., 1., -1., 1., -1., 1., -1., 1.])
        self.noise = self.h * self.s * self.c
        self.design = np.column_stack((np.ones(8), self.h, self.s, self.c))
        self.assertEqual(np.linalg.matrix_rank(self.design), 4)

    def test_response_scale_slope_is_zero_for_constant_response(self):
        # This is an OLS coefficient on a different scale, not imputation of
        # the old undefined standardized coefficient.
        for constant in (0.0, 0.2):
            y = np.full(8, constant)
            fit = np.linalg.lstsq(self.design, y, rcond=None)[0]
            np.testing.assert_allclose(fit[1:], 0.0, atol=1e-14)
            self.assertEqual(float(np.std(y)), 0.0)

    def test_exact_relation_to_original_standardized_coefficient(self):
        y = (0.012 + 0.003*self.h - 0.002*self.s
             + 0.001*self.c + 0.0005*self.noise)
        sd = float(np.std(y, ddof=0))
        original = np.linalg.lstsq(self.design, (y-y.mean())/sd, rcond=None)[0][1]
        response_scale = np.linalg.lstsq(self.design, y, rcond=None)[0][1]
        self.assertAlmostEqual(response_scale, 0.003, places=14)
        self.assertAlmostEqual(response_scale, sd*original, places=14)

    def test_partial_projection_and_full_ols_agree(self):
        # r is the part of standardized host history left after controls.
        host = self.h + 0.3*self.s
        host = (host-host.mean())/host.std(ddof=0)
        controls = np.column_stack((np.ones(8), self.s, self.c))
        residual_host = host-controls@np.linalg.lstsq(controls, host, rcond=None)[0]
        y = 0.4 + 0.005*host + 0.002*self.s - 0.003*self.c + 0.0002*self.noise
        projected = float(residual_host@y/(residual_host@host))
        full = np.linalg.lstsq(np.column_stack((controls[:,0], host, self.s, self.c)),
                               y, rcond=None)[0][1]
        self.assertAlmostEqual(projected, full, places=14)

    def test_standardized_slope_has_no_unique_extension_at_zero_response(self):
        # Approaching y=0 from opposite directions gives opposite standardized
        # slopes, while both response-scale slopes converge to zero.
        for eps in (1.0, 1e-3, 1e-9):
            plus = eps*self.h
            minus = -eps*self.h
            bp = np.linalg.lstsq(self.design, plus/plus.std(), rcond=None)[0][1]
            bm = np.linalg.lstsq(self.design, minus/minus.std(), rcond=None)[0][1]
            gp = np.linalg.lstsq(self.design, plus, rcond=None)[0][1]
            gm = np.linalg.lstsq(self.design, minus, rcond=None)[0][1]
            self.assertAlmostEqual(bp, 1.0, places=13)
            self.assertAlmostEqual(bm, -1.0, places=13)
            self.assertAlmostEqual(gp, eps, places=13)
            self.assertAlmostEqual(gm, -eps, places=13)
        with np.errstate(invalid='ignore', divide='ignore'):
            standardized_zero = np.zeros(8)/np.std(np.zeros(8))
        self.assertTrue(np.all(np.isnan(standardized_zero)))

    def test_average_response_scale_slope_is_a_different_estimand(self):
        # Illustrative synthetic species, not any observed TTF coefficients.
        standardized = np.array([1.0, -1.0])
        response_sd = np.array([0.01, 0.001])
        self.assertEqual(float(standardized.mean()), 0.0)
        response_scale = standardized*response_sd
        self.assertAlmostEqual(float(response_scale.mean()), 0.0045, places=14)
        all_four = np.r_[response_scale, 0.0, 0.0]
        self.assertAlmostEqual(float(all_four.mean()), 0.00225, places=14)

    def test_original_confirmatory_terminal_result_is_preserved(self):
        root = Path(__file__).resolve().parents[1]
        path = root/'benchmarks/frozen/historical_host_memory_259_one_shot_empirical_result_v0.1.json'
        raw = path.read_bytes()
        self.assertEqual(hashlib.sha256(raw).hexdigest(),
                         '5d994a76f7f20a778f0394f545c08a729182ce32298503c128d95c00da29f965')
        receipt = json.loads(raw)
        self.assertEqual(receipt['status'], 'NOT_EVALUABLE_EMPIRICAL_RESPONSE')
        self.assertEqual(receipt['n_species_frozen'], 259)
        self.assertEqual(receipt['undefined_species'], 4)
        self.assertIsNone(receipt['primary_statistic'])
        self.assertIsNone(receipt['monte_carlo_p'])
        self.assertIs(receipt['empirical_identity_opened'], True)


if __name__ == '__main__':
    unittest.main()
