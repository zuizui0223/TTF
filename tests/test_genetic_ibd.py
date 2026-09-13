import numpy as np
import pytest

from ttf.core import knn_edges, spearman_rho
from ttf.genetic_ibd import (
    crossfit_ibd_residuals,
    crossfit_ibd_residuals_prepared,
    prepare_crossfit_ibd_design,
)


def _random_edge_geometry(seed: int = 20260913):
    rng = np.random.default_rng(seed)
    coordinates = rng.normal(size=(20, 2))
    edges = knn_edges(coordinates, k=4)
    geographic = np.linalg.norm(
        coordinates[edges[:, 1]] - coordinates[edges[:, 0]],
        axis=1,
    )
    return rng, coordinates, edges, geographic


def test_crossfit_ibd_removes_noiseless_strictly_monotone_ibd_exactly():
    _, _, edges, geographic = _random_edge_geometry()
    genetic = np.exp(geographic)
    result = crossfit_ibd_residuals(genetic, geographic, edges)
    assert np.array_equal(result.residual, np.zeros_like(result.residual))
    assert np.array_equal(result.residual_turnover, np.full_like(result.residual, 0.5))


def test_crossfit_ibd_removes_tie_rich_monotone_ibd_exactly():
    _, _, edges, geographic = _random_edge_geometry(17)
    geographic = np.round(geographic, 1)
    genetic = 0.137 * geographic
    result = crossfit_ibd_residuals(genetic, geographic, edges)
    assert np.array_equal(result.residual, np.zeros_like(result.residual))
    assert np.array_equal(result.residual_turnover, np.full_like(result.residual, 0.5))
    assert np.array_equal(result.observed_rank_fraction, result.expected_rank_fraction)


def test_crossfit_ibd_near_tied_direct_scalar_relation_is_exact_zero():
    _, _, edges, geographic = _random_edge_geometry(41)
    # Force two distances to be adjacent floating-point neighbours. The direct
    # positive scalar relation must retain the same weak order and be evaluated
    # as exactly pure rank-IBD rather than amplifying roundoff through rank01.
    geographic = geographic.copy()
    geographic[1] = np.nextafter(geographic[0], np.inf)
    genetic = geographic.copy()
    result = crossfit_ibd_residuals(genetic, geographic, edges)
    assert np.array_equal(result.residual, np.zeros_like(result.residual))


def test_prepared_and_wrapper_crossfit_are_exactly_equivalent_nonmonotone():
    rng, _, edges, geographic = _random_edge_geometry(55)
    genetic = geographic + rng.normal(0.0, 0.07, len(geographic))
    genetic -= float(genetic.min()) - 0.01
    direct = crossfit_ibd_residuals(genetic, geographic, edges)
    design = prepare_crossfit_ibd_design(geographic, edges)
    cached = crossfit_ibd_residuals_prepared(genetic, design)
    assert np.array_equal(direct.n_training_edges, cached.n_training_edges)
    assert np.allclose(direct.residual, cached.residual, atol=0.0, rtol=0.0)
    assert np.allclose(direct.residual_turnover, cached.residual_turnover, atol=0.0, rtol=0.0)
    assert np.allclose(direct.expected_rank_fraction, cached.expected_rank_fraction, atol=0.0, rtol=0.0)


def test_crossfit_ibd_is_invariant_to_strictly_monotone_rescaling():
    rng, _, edges, geographic = _random_edge_geometry(33)
    genetic = 0.4 * geographic + rng.normal(0.0, 0.03, len(geographic))
    genetic -= float(genetic.min()) - 0.01

    original = crossfit_ibd_residuals(genetic, geographic, edges)
    transformed = crossfit_ibd_residuals(
        np.exp(genetic),
        geographic**3 + 7.0,
        edges,
    )
    assert np.allclose(original.residual, transformed.residual, atol=1e-12, rtol=0.0)
    assert np.allclose(
        original.residual_turnover,
        transformed.residual_turnover,
        atol=1e-12,
        rtol=0.0,
    )


def test_target_edge_nuisance_fit_ignores_outcomes_touching_its_localities():
    rng, _, edges, geographic = _random_edge_geometry(77)
    genetic = geographic + rng.normal(0.0, 0.05, len(geographic))
    genetic -= float(genetic.min()) - 0.01
    target = 0
    left, right = edges[target]

    altered = genetic.copy()
    incident = (
        (edges[:, 0] == left)
        | (edges[:, 1] == left)
        | (edges[:, 0] == right)
        | (edges[:, 1] == right)
    )
    incident[target] = False
    altered[incident] += 1000.0

    baseline = crossfit_ibd_residuals(genetic, geographic, edges)
    perturbed = crossfit_ibd_residuals(altered, geographic, edges)
    assert np.isclose(
        baseline.residual[target],
        perturbed.residual[target],
        atol=1e-12,
        rtol=0.0,
    )
    assert np.isclose(
        baseline.expected_rank_fraction[target],
        perturbed.expected_rank_fraction[target],
        atol=1e-12,
        rtol=0.0,
    )


def test_crossfit_ibd_preserves_barrier_departure_from_distance_expectation():
    rng, coordinates, edges, geographic = _random_edge_geometry(91)
    barrier = (
        coordinates[edges[:, 0], 0] * coordinates[edges[:, 1], 0] < 0.0
    ).astype(float)
    genetic = geographic + 1.2 * barrier + rng.normal(0.0, 0.02, len(geographic))
    genetic -= float(genetic.min()) - 0.01

    result = crossfit_ibd_residuals(genetic, geographic, edges)
    assert spearman_rho(result.residual_turnover, barrier) > 0.6
    assert float(result.residual[barrier == 1].mean()) > float(
        result.residual[barrier == 0].mean()
    )


def test_crossfit_ibd_fails_when_endpoint_disjoint_training_set_is_too_small():
    edges = np.asarray([[0, 1], [0, 2], [1, 2]], dtype=np.int64)
    geographic = np.asarray([1.0, 1.2, 0.8])
    genetic = np.asarray([0.1, 0.2, 0.15])
    with pytest.raises(ValueError, match="not estimable"):
        crossfit_ibd_residuals(genetic, geographic, edges, min_training_edges=3)
