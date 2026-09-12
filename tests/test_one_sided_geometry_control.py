import numpy as np

from ttf.core import average_ranks
from ttf.one_sided_geometry_control import (
    centered_rank,
    rank_residual,
    residual_rank_correlation,
)
from ttf.private_geometry_control import partial_spearman_controls


def test_rank_residual_is_orthogonal_to_nuisance_ranks():
    z1 = np.linspace(-2, 2, 101)
    z2 = np.sin(np.linspace(0, 4 * np.pi, 101))
    y = 3 * z1 - 2 * z2 + 0.1 * np.cos(z1)
    r = rank_residual(y, [z1, z2])
    assert abs(float(r.mean())) < 1e-12
    for z in (z1, z2):
        rz = average_ranks(z)
        rz -= rz.mean()
        assert abs(float(np.dot(r, rz))) < 1e-9


def test_no_nuisance_matches_spearman_rank_correlation():
    rng = np.random.default_rng(33)
    x = rng.normal(size=120)
    y = 0.3 * x + rng.normal(size=120)
    got = residual_rank_correlation(x, y)
    rx = centered_rank(x)
    ry = centered_rank(y)
    expected = np.dot(rx, ry) / np.sqrt(np.dot(rx, rx) * np.dot(ry, ry))
    assert np.isclose(got, expected, atol=1e-12, rtol=0)


def test_two_sided_controls_match_existing_partial_spearman():
    rng = np.random.default_rng(44)
    x = rng.normal(size=150)
    y = rng.normal(size=150)
    z1 = rng.normal(size=150)
    z2 = rng.normal(size=150)
    old = partial_spearman_controls(x, y, [z1, z2])
    new = residual_rank_correlation(
        x, y,
        prediction_nuisances=[z1, z2],
        target_nuisances=[z1, z2],
    )
    assert np.isclose(old, new, atol=1e-12, rtol=0)


def test_prediction_only_differs_from_two_sided_when_target_contains_signal_aligned_with_nuisance():
    z = np.linspace(-3, 3, 201)
    signal = np.sin(z)
    prediction = signal + 0.5 * z
    target = signal + 0.5 * z
    one = residual_rank_correlation(
        prediction, target, prediction_nuisances=[z], target_nuisances=[]
    )
    two = residual_rank_correlation(
        prediction, target, prediction_nuisances=[z], target_nuisances=[z]
    )
    assert np.isfinite(one) and np.isfinite(two)
    assert not np.isclose(one, two)
