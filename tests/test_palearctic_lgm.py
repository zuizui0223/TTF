import numpy as np

from ttf.palearctic_lgm import (
    climate_features,
    fit_climate_scaling,
    fit_weighted_ridge_logistic,
    predict_relative_suitability,
    relation_matrix,
    schoener_d,
    target_group_background_cells,
)


def test_target_group_background_is_order_invariant_and_deduplicated():
    x=np.asarray([
        [10.01,50.01],
        [10.02,50.02],
        [11.01,51.01],
        [12.01,52.01],
    ])
    a=target_group_background_cells(x,resolution_degrees=.25)
    b=target_group_background_cells(x[::-1],resolution_degrees=.25)
    np.testing.assert_allclose(a,b)
    assert len(a)==3


def test_climate_scaling_features_are_finite_and_clamped():
    rng=np.random.default_rng(2)
    background=rng.normal(size=(500,4))
    scaling=fit_climate_scaling(background)
    features=climate_features(
        np.vstack((background[:5],np.full((1,4),1e6))),
        scaling,
    )
    assert features.shape==(6,8)
    assert np.isfinite(features).all()
    assert np.max(np.abs(features[-1,:4])) < 10


def test_ridge_logistic_recovers_positive_climate_direction():
    rng=np.random.default_rng(4)
    background=rng.normal(size=(1000,4))
    presence=rng.normal(size=(100,4))
    presence[:,0]+=2.0
    scaling=fit_climate_scaling(background)
    p=climate_features(presence,scaling)
    b=climate_features(background,scaling)
    fit=fit_weighted_ridge_logistic(p,b)
    assert fit.converged
    assert fit.coefficients[1] > 0


def test_suitability_normalization_and_overlap():
    rng=np.random.default_rng(7)
    background=rng.normal(size=(1000,4))
    p1=rng.normal(size=(80,4)); p1[:,0]+=1.0
    p2=rng.normal(size=(80,4)); p2[:,0]-=1.0
    scaling=fit_climate_scaling(background)
    bg=climate_features(background,scaling)
    f1=fit_weighted_ridge_logistic(climate_features(p1,scaling),bg)
    f2=fit_weighted_ridge_logistic(climate_features(p2,scaling),bg)
    grid=climate_features(rng.normal(size=(200,4)),scaling)
    s1=predict_relative_suitability(grid,f1)
    s2=predict_relative_suitability(grid,f2)
    assert np.isclose(s1.sum(),1.0)
    assert np.isclose(s2.sum(),1.0)
    assert np.isclose(schoener_d(s1,s1),1.0)
    assert 0 <= schoener_d(s1,s2) <= 1
    mat=relation_matrix(np.vstack((s1,s2)))
    assert mat.shape==(2,2)
    assert np.allclose(mat,mat.T)
    assert np.allclose(np.diag(mat),1.0)
