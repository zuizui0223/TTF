from __future__ import annotations

import numpy as np

from ttf.rank_factorization import standardized_rank_direction


def test_rank_direction_is_odd_under_global_sign_reversal_with_ties() -> None:
    values = np.array([4.0, 4.0, -1.0, 2.0, -3.0, -3.0])
    assert np.allclose(
        standardized_rank_direction(-values),
        -standardized_rank_direction(values),
        atol=1e-15,
        rtol=0.0,
    )


def test_centrally_symmetric_support_has_zero_expected_rank_direction() -> None:
    a = np.array([8.0, 2.0, -1.0, -5.0, 3.0])
    b = np.array([1.0, 7.0, -4.0, 2.0, -9.0])
    support = np.vstack([a, -a, b, -b])

    # The coordinates have strongly unequal variances, so this is not an
    # equal-variance construction.  Central symmetry alone forces cancellation.
    coordinate_variances = np.var(support, axis=0)
    assert float(np.max(coordinate_variances) - np.min(coordinate_variances)) > 10.0

    mean_direction = np.mean(
        np.vstack([standardized_rank_direction(row) for row in support]),
        axis=0,
    )
    assert np.allclose(mean_direction, 0.0, atol=1e-15, rtol=0.0)


def test_zero_raw_mean_does_not_force_zero_expected_rank_direction() -> None:
    support = np.array(
        [
            [-3.0, 3.0, 1.0, 4.0],
            [0.0, 1.0, 4.0, 2.0],
            [3.0, -4.0, -5.0, -6.0],
        ]
    )
    assert np.allclose(np.mean(support, axis=0), 0.0, atol=0.0, rtol=0.0)

    mean_direction = np.mean(
        np.vstack([standardized_rank_direction(row) for row in support]),
        axis=0,
    )
    assert np.linalg.norm(mean_direction) > 0.25


def test_constant_vector_obeys_odd_identity_at_zero() -> None:
    values = np.full(5, 7.0)
    assert np.array_equal(standardized_rank_direction(values), np.zeros(5))
    assert np.array_equal(standardized_rank_direction(-values), np.zeros(5))
