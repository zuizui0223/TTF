import numpy as np

from ttf import (
    centered_species_bootstrap_mean_test,
    heldout_species_bootstrap_test,
    simulate_circular_boundary_world,
    split_species,
)


def test_centered_species_bootstrap_is_deterministic_and_null_centered():
    scores = np.array([-0.30, -0.20, -0.10, 0.10, 0.20, 0.30])
    a = centered_species_bootstrap_mean_test(scores, n_bootstrap=999, seed=17)
    b = centered_species_bootstrap_mean_test(scores, n_bootstrap=999, seed=17)
    assert np.array_equal(a.null_means, b.null_means)
    assert np.array_equal(a.null_studentized, b.null_studentized)
    assert a.p_value == b.p_value
    assert abs(float(a.null_means.mean())) < 0.01
    assert a.p_value >= 0.25


def test_centered_species_bootstrap_detects_positive_species_mean():
    scores = np.array([0.10, 0.18, 0.24, 0.31, 0.42, 0.55, 0.61, 0.73])
    result = centered_species_bootstrap_mean_test(scores, n_bootstrap=1999, seed=91)
    assert result.observed_mean > 0.0
    assert result.observed_studentized > 3.0
    assert result.p_value <= 0.01


def test_heldout_species_bootstrap_matches_primary_transfer_mean():
    world = simulate_circular_boundary_world(
        n_species=16,
        records_per_species=28,
        shared_fraction=1.0,
        amplitude=3.0,
        noise_sd=0.35,
        transition_width=0.15,
        seed=123,
    )
    train, evaluation = split_species(
        [sample.species for sample in world.samples],
        eval_fraction=0.5,
        seed=456,
    )
    result = heldout_species_bootstrap_test(
        world.samples,
        train_species=train,
        eval_species=evaluation,
        k=4,
        bandwidth=0.2,
        n_bootstrap=199,
        seed=789,
    )
    assert result.observed.n_eval_species == len(evaluation)
    assert np.isclose(
        result.observed.statistic,
        np.mean(result.species_scores),
        atol=1e-12,
        rtol=0.0,
    )
    assert abs(float(result.null_statistics.mean())) < 0.02
    assert result.observed.statistic > 0.05
    assert result.p_value <= 0.05


def test_heldout_species_inference_does_not_modify_input_samples():
    world = simulate_circular_boundary_world(
        n_species=12,
        records_per_species=20,
        shared_fraction=0.0,
        amplitude=2.0,
        seed=999,
    )
    before = {
        sample.species: (sample.coordinates.copy(), sample.trait.copy())
        for sample in world.samples
    }
    train, evaluation = split_species(
        [sample.species for sample in world.samples], eval_fraction=0.5, seed=12
    )
    heldout_species_bootstrap_test(
        world.samples,
        train_species=train,
        eval_species=evaluation,
        bandwidth=0.2,
        n_bootstrap=99,
        seed=13,
    )
    for sample in world.samples:
        coords, trait = before[sample.species]
        assert np.array_equal(sample.coordinates, coords)
        assert np.array_equal(sample.trait, trait)
