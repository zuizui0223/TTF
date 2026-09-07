import numpy as np

from ttf.core import SpeciesEdges, average_ranks
from ttf.orthogonal import orthogonalize_rank_predictor, semi_partial_rank_score


def test_rank_predictor_is_orthogonal_to_active_geometry_covariates():
    opportunity = np.linspace(0.1, 2.0, 30)
    edge_length = np.linspace(2.0, 0.2, 30) + 0.1 * np.sin(np.arange(30))
    predictor = 4.0 * opportunity + 2.0 * edge_length + 0.2 * np.cos(np.arange(30))
    covariates = np.column_stack([opportunity, edge_length])
    out = orthogonalize_rank_predictor(predictor, covariates)
    assert out.n_active_covariates == 2
    assert out.nuisance_rank_r2 > 0.9
    assert abs(float(out.residual.mean())) < 1e-12
    for j in range(covariates.shape[1]):
        ranked = average_ranks(covariates[:, j])
        centered = ranked - ranked.mean()
        assert abs(float(np.dot(out.residual, centered))) < 1e-9


def test_semi_partial_score_removes_geometry_only_association():
    opportunity = np.linspace(0.0, 1.0, 40)
    predictor = opportunity + 0.02 * np.sin(np.arange(40))
    target = opportunity + 0.02 * np.cos(np.arange(40))
    score, _ = semi_partial_rank_score(predictor, target, opportunity)
    assert abs(score) < 0.25


def test_semi_partial_score_retains_component_not_explained_by_geometry():
    x = np.linspace(-1.0, 1.0, 50)
    opportunity = x * x
    shared = np.sin(4.0 * x)
    predictor = 3.0 * opportunity + shared
    target = 2.0 * opportunity + shared
    score, _ = semi_partial_rank_score(predictor, target, opportunity)
    assert score > 0.45
