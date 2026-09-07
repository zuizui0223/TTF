import numpy as np

from ttf.controls import paired_opportunity_control_test, score_paired_opportunity_control
from ttf.core import build_species_edges, split_species
from ttf.simulate import simulate_circular_boundary_world
from ttf.transfer import prepare_transfer


def test_opportunity_control_is_trait_free_given_geometry():
    world = simulate_circular_boundary_world(
        n_species=8,
        records_per_species=24,
        shared_fraction=1.0,
        amplitude=2.0,
        seed=101,
    )
    edges = {sample.species: build_species_edges(sample, k=4) for sample in world.samples}
    train, evaluation = split_species(tuple(edges), eval_fraction=0.75, seed=3)
    prepared = prepare_transfer(
        [edges[name] for name in train],
        [edges[name] for name in evaluation],
        bandwidth=0.2,
    )
    original_support = {name: prepared.eval_opportunity[name].copy() for name in evaluation}
    altered = {
        name: np.linspace(0.01, 0.99, edges[name].n_edges)
        for name in train
    }
    score_paired_opportunity_control(
        prepared,
        train_turnover=altered,
        eval_turnover={name: edges[name].turnover for name in evaluation},
    )
    for name in evaluation:
        assert np.array_equal(prepared.eval_opportunity[name], original_support[name])


def test_paired_control_statistic_equals_boundary_minus_opportunity_macroaverage():
    world = simulate_circular_boundary_world(
        n_species=8,
        records_per_species=28,
        shared_fraction=1.0,
        amplitude=2.5,
        noise_sd=0.5,
        seed=103,
    )
    edges = {sample.species: build_species_edges(sample, k=4) for sample in world.samples}
    train, evaluation = split_species(tuple(edges), eval_fraction=0.75, seed=5)
    result = paired_opportunity_control_test(
        [edges[name] for name in train],
        [edges[name] for name in evaluation],
        bandwidth=0.2,
        n_bootstrap=199,
        seed=7,
    )
    deltas = []
    for name, value in result.observed.species_scores.items():
        assert np.isclose(
            value,
            result.observed.boundary_scores[name] - result.observed.opportunity_scores[name],
            atol=1e-12,
        )
        if np.isfinite(value):
            deltas.append(value)
    assert len(deltas) >= 6
    assert np.isclose(result.observed.statistic, np.mean(deltas), atol=1e-12)
    assert 0.0 < result.p_value <= 1.0
