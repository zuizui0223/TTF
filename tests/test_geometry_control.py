import numpy as np

from ttf.geometry_control import (
    length_orthogonalized_turnover,
    partial_spearman_rho,
    spearman_length_association,
)


def test_length_orthogonalized_turnover_removes_perfect_monotone_length_signal():
    length = np.arange(1.0, 21.0)
    turnover = np.linspace(0.01, 0.99, len(length))
    residual = length_orthogonalized_turnover(turnover, length)
    assert np.allclose(residual, 0.0, atol=1e-12, rtol=0.0)


def test_length_orthogonalization_is_invariant_to_monotone_length_rescaling():
    rng = np.random.default_rng(20260908)
    length = np.sort(rng.uniform(0.1, 5.0, 30))
    turnover = rng.normal(size=30)
    a = length_orthogonalized_turnover(turnover, length)
    b = length_orthogonalized_turnover(turnover, np.exp(length))
    assert np.allclose(a, b, atol=1e-12, rtol=0.0)


def test_length_orthogonalized_turnover_has_zero_linear_rank_covariance():
    rng = np.random.default_rng(7)
    length = rng.lognormal(size=40)
    turnover = 0.7 * np.argsort(np.argsort(length)) + rng.normal(0, 5, 40)
    residual = length_orthogonalized_turnover(turnover, length)
    length_rank = np.argsort(np.argsort(length)).astype(float)
    length_rank -= length_rank.mean()
    assert abs(float(np.dot(residual, length_rank))) < 1e-10
    assert abs(float(residual.mean())) < 1e-12


def test_partial_spearman_removes_common_edge_length_path():
    rng = np.random.default_rng(91)
    length = np.linspace(0.0, 1.0, 80)
    prediction = length + rng.normal(0.0, 0.03, len(length))
    turnover = length + rng.normal(0.0, 0.03, len(length))
    raw_proxy = spearman_length_association(turnover, length)
    partial = partial_spearman_rho(prediction, turnover, length)
    assert raw_proxy > 0.9
    assert abs(partial) < 0.25


def test_partial_spearman_preserves_shared_signal_not_explained_by_length():
    rng = np.random.default_rng(19)
    length = rng.uniform(size=120)
    shared = rng.normal(size=120)
    prediction = shared + 0.8 * length + rng.normal(0.0, 0.1, 120)
    turnover = shared + 0.8 * length + rng.normal(0.0, 0.1, 120)
    partial = partial_spearman_rho(prediction, turnover, length)
    assert partial > 0.8
