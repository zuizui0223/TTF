import numpy as np
import pytest

from ttf.historical_host_memory import (
    equal_species_mean,
    historical_analog_memory,
    nearest_cloud_distance,
    partial_host_memory_beta,
    residualized_unique_fraction,
    standardize_from_reference,
)


def test_standardization_uses_reference_only():
    ref = np.array([[0.0, 0.0], [2.0, 4.0], [4.0, 8.0]])
    a, b = standardize_from_reference(
        ref, np.array([[2.0, 4.0]]), np.array([[4.0, 8.0]])
    )
    assert np.allclose(a, [[0.0, 0.0]])
    assert np.all(b > 0)


def test_historical_memory_is_positive_when_lgm_is_farther_from_host_climate():
    host = np.array([[0.0, 0.0]])
    current = np.zeros((2, 5, 2))
    lgm = np.array([[[2.0, 0.0]] * 5, [[4.0, 0.0]] * 5])
    memory, d_current, d_lgm = historical_analog_memory(current, lgm, host)
    assert np.allclose(d_current, [0.0, 0.0])
    assert np.allclose(d_lgm, [2.0, 4.0])
    assert np.allclose(memory, [2.0, 4.0])


def test_nearest_cloud_uses_nearest_analog_not_centroid():
    cloud = np.array([[0.0, 0.0], [10.0, 0.0]])
    points = np.array([[1.0, 0.0], [9.0, 0.0], [5.0, 0.0]])
    assert np.allclose(nearest_cloud_distance(points, cloud), [1.0, 1.0, 5.0])


def test_partial_beta_recovers_host_memory_after_controls():
    rng = np.random.default_rng(4)
    n = 300
    host = rng.normal(size=n)
    self_memory = 0.45 * host + rng.normal(size=n)
    current = -0.35 * host + rng.normal(size=n)
    turnover = (
        0.8 * host
        + 0.5 * self_memory
        - 0.4 * current
        + rng.normal(scale=0.25, size=n)
    )
    beta = partial_host_memory_beta(turnover, host, self_memory, current)
    assert beta > 0.5


def test_unique_fraction_drops_when_primary_is_explained_by_controls():
    x = np.linspace(-2, 2, 100)
    y = x + 0.01 * np.sin(np.arange(100))
    fraction = residualized_unique_fraction(y, np.column_stack([x]))
    assert fraction < 0.01


def test_equal_species_mean_does_not_weight_dense_species():
    assert equal_species_mean([1.0, -1.0, 0.5]) == pytest.approx(1.0 / 6.0)