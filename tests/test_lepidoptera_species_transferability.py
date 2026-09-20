from __future__ import annotations

import numpy as np

from ttf.lepidoptera_species_transferability import (
    crossvalidated_species_property_gain,
    decompose_species_transferability,
    deterministic_pair_folds,
)


def synthetic_pairs(seed: int = 4):
    rng=np.random.default_rng(seed)
    sources=[f"s{i:02d}" for i in range(18)]
    targets=[f"t{i:02d}" for i in range(14)]
    source=[]; target=[]; geometry=[]; y=[]
    u=dict(zip(sources,rng.normal(0,.35,len(sources))))
    v=dict(zip(targets,rng.normal(0,.22,len(targets))))
    beta=np.asarray([.12,-.08,.05,.03,-.06])
    for s in sources:
        for t in targets:
            g=rng.normal(size=5)
            source.append(s); target.append(t); geometry.append(g)
            y.append(.1+g@beta+u[s]+v[t]+rng.normal(0,.12))
    return np.asarray(y),np.asarray(geometry),source,target


def test_decomposition_recovers_nonzero_source_and_target_components():
    y,g,s,t=synthetic_pairs()
    result=decompose_species_transferability(y,g,s,t)
    assert result.source_variance>0.03
    assert result.target_variance>0.01
    assert result.residual_variance>0
    assert np.isclose(
        result.source_fraction+result.target_fraction+result.residual_fraction,
        1.0,
    )
    assert len(result.source_effect)==18
    assert len(result.target_effect)==14


def test_species_property_improves_pair_prediction_when_species_effects_exist():
    y,g,s,t=synthetic_pairs()
    result=crossvalidated_species_property_gain(y,g,s,t,folds=10)
    assert result.species_property_rmse<result.geometry_only_rmse
    assert result.rmse_improvement>0
    assert result.species_property_correlation>result.geometry_only_correlation


def test_pair_folds_are_deterministic_and_pair_label_based():
    y,g,s,t=synthetic_pairs()
    a=deterministic_pair_folds(s,t,folds=10)
    b=deterministic_pair_folds(list(s),list(t),folds=10)
    assert np.array_equal(a,b)
    assert set(a)<=set(range(10))
