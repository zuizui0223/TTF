import numpy as np

from ttf.core import knn_edges, rank01
from ttf.geometry import SpeciesGeometry
from ttf.private_geometry_control import (
    partial_spearman_controls,
    pooled_affine_frame,
    random_private_transition_propensity,
    rank_residualized_turnover,
)


def test_random_private_propensity_is_deterministic_and_geometry_only():
    x = np.column_stack([np.linspace(-2, 2, 21), np.sin(np.linspace(-2, 2, 21))])
    g = SpeciesGeometry("sp", x)
    center, scale, active = pooled_affine_frame([g])
    edges = knn_edges(x, k=3)
    a = random_private_transition_propensity(
        x, edges, center=center, scale=scale, active=active,
        n_directions=256, seed=17,
    )
    b = random_private_transition_propensity(
        x, edges, center=center, scale=scale, active=active,
        n_directions=256, seed=17,
    )
    assert a.shape == (len(edges),)
    assert np.isfinite(a).all()
    assert np.all(a >= 0)
    assert np.allclose(a, b, atol=0, rtol=0)


def test_rank_residualization_removes_frozen_nuisance_span():
    n = 80
    z1 = np.linspace(-1, 1, n)
    z2 = np.cos(np.linspace(0, 3 * np.pi, n))
    y = 4 * rank01(z1) - 2 * rank01(z2) + np.linspace(0, 0.01, n)
    residual = rank_residualized_turnover(y, [z1, z2])
    assert abs(float(residual.mean())) < 1e-12
    for z in (z1, z2):
        rz = rank01(z) - float(rank01(z).mean())
        assert abs(float(np.dot(residual, rz))) < 1e-10


def test_partial_spearman_controls_removes_shared_geometry_nuisance():
    z = np.linspace(-3, 3, 101)
    x = z + 0.03 * np.sin(5 * z)
    y = z + 0.03 * np.cos(4 * z)
    raw = np.corrcoef(rank01(x), rank01(y))[0, 1]
    partial = partial_spearman_controls(x, y, [z])
    assert raw > 0.95
    assert abs(partial) < 0.2
