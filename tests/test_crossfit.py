import numpy as np

from ttf.core import build_species_edges
from ttf.crossfit import balanced_species_folds, crossfit_species_bootstrap_test
from ttf.simulate import simulate_circular_boundary_world


def test_balanced_species_folds_partition_every_species_once():
    labels = [f"sp{i}" for i in range(9)]
    folds = balanced_species_folds(labels, n_folds=3, seed=17)
    assert [len(fold) for fold in folds] == [3, 3, 3]
    flat = [name for fold in folds for name in fold]
    assert len(flat) == len(set(flat)) == 9
    assert set(flat) == set(labels)
    assert folds == balanced_species_folds(labels, n_folds=3, seed=17)


def test_crossfit_scores_every_species_once_and_bootstraps_macro_mean():
    world = simulate_circular_boundary_world(
        n_species=9,
        records_per_species=24,
        shared_fraction=1.0,
        amplitude=2.5,
        noise_sd=0.4,
        seed=29,
    )
    edge_sets = tuple(build_species_edges(sample, k=4) for sample in world.samples)
    result, folds = crossfit_species_bootstrap_test(
        edge_sets,
        n_folds=3,
        fold_seed=31,
        bandwidth=0.2,
        n_bootstrap=199,
        seed=37,
    )
    assert len(folds) == 3
    assert set(result.observed.species_scores) == {sample.species for sample in world.samples}
    assert result.observed.n_eval_species == 9
    assert np.isclose(result.observed.statistic, np.mean(result.species_scores), atol=1e-12)
    assert 0.0 < result.p_value <= 1.0
