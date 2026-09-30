import numpy as np

from ttf.lgm_climate_envelope import (
    climatic_relation,
    current_lgm_suitability,
    fit_common_climate_space,
    gaussian_kde_density,
    normalized_suitability,
    scott_factor,
)


def synthetic_climate(seed=3):
    rng=np.random.default_rng(seed)
    pooled=rng.normal(size=(200,4))
    species=pooled[:40]+np.array([0.3,-0.2,0.1,0.0])
    current=rng.normal(size=(120,4))
    lgm=current+np.array([1.0,-0.6,0.4,0.2])
    return pooled,species,current,lgm


def test_four_dimensional_scott_factor_matches_frozen_contract():
    assert np.isclose(scott_factor(100,4),100**(-1/8))


def test_suitability_is_normalized_without_thresholding():
    pooled,species,current,_=synthetic_climate()
    space=fit_common_climate_space(pooled,axes=4)
    suitability=normalized_suitability(species,current,space)
    assert suitability.shape==(len(current),)
    assert np.isfinite(suitability).all()
    assert np.all(suitability>=0)
    assert np.isclose(suitability.sum(),1.0)


def test_identical_climate_projection_has_identical_relation():
    pooled,species,current,_=synthetic_climate()
    space=fit_common_climate_space(pooled,axes=4)
    now,lgm=current_lgm_suitability(species,current,current.copy(),space)
    assert np.allclose(now,lgm)
    assert np.isclose(climatic_relation(now,lgm),1.0)


def test_lgm_shift_changes_climatic_suitability_surface():
    pooled,species,current,lgm=synthetic_climate()
    space=fit_common_climate_space(pooled,axes=4)
    now,past=current_lgm_suitability(species,current,lgm,space)
    assert climatic_relation(now,past)<0.95


def test_kde_rejects_dimension_drift():
    x=np.zeros((4,4))
    q=np.zeros((3,3))
    try:
        gaussian_kde_density(x,q)
    except ValueError:
        pass
    else:
        raise AssertionError("dimension drift must fail closed")
