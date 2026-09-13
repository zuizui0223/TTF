from __future__ import annotations

import numpy as np

from ttf.core import spearman_rho
from ttf.rank_factorization import (
    independent_rank_factorization,
    spearman_rank_inner_product,
    standardized_rank_direction,
)


def test_rank_inner_product_matches_core_spearman_with_ties() -> None:
    x = np.array([3.0, 1.0, 1.0, 7.0, 5.0, 5.0])
    y = np.array([0.0, 2.0, 1.0, 1.0, 4.0, 3.0])
    assert np.isclose(
        spearman_rank_inner_product(x, y),
        spearman_rho(x, y),
        atol=1e-15,
        rtol=0.0,
    )


def test_constant_direction_matches_zero_score_contract() -> None:
    constant = np.ones(5)
    varying = np.arange(5.0)
    assert np.array_equal(standardized_rank_direction(constant), np.zeros(5))
    assert spearman_rank_inner_product(constant, varying) == 0.0
    assert spearman_rank_inner_product(varying, constant) == 0.0


def test_independent_finite_support_expectation_factorizes_exactly() -> None:
    left = np.array(
        [
            [3.0, 1.0, 2.0, 5.0, 4.0],
            [1.0, 1.0, 4.0, 3.0, 2.0],
            [4.0, 2.0, 2.0, 1.0, 5.0],
        ]
    )
    right = np.array(
        [
            [0.0, 1.0, 0.0, 1.0, 1.0],
            [1.0, 0.0, 1.0, 0.0, 0.0],
            [2.0, 1.0, 4.0, 3.0, 0.0],
            [7.0, 7.0, 7.0, 7.0, 7.0],
        ]
    )
    lw = np.array([0.2, 0.3, 0.5])
    rw = np.array([0.1, 0.25, 0.4, 0.25])

    factorized = independent_rank_factorization(
        left,
        right,
        left_weights=lw,
        right_weights=rw,
    )
    direct = 0.0
    for i in range(len(left)):
        for j in range(len(right)):
            direct += lw[i] * rw[j] * spearman_rho(left[i], right[j])

    assert np.isclose(factorized.expected_spearman, direct, atol=1e-15, rtol=0.0)


def test_factorization_is_zero_when_one_expected_rank_direction_is_zero() -> None:
    left = np.array(
        [
            [0.0, 1.0, 2.0, 3.0],
            [3.0, 2.0, 1.0, 0.0],
        ]
    )
    right = np.array(
        [
            [0.0, 0.0, 1.0, 1.0],
            [1.0, 1.0, 0.0, 0.0],
        ]
    )
    result = independent_rank_factorization(left, right)
    assert np.allclose(result.mean_left_direction, 0.0, atol=1e-15, rtol=0.0)
    assert np.allclose(result.mean_right_direction, 0.0, atol=1e-15, rtol=0.0)
    assert abs(result.expected_spearman) <= 1e-15
