import numpy as np

from ttf import (
    build_species_edges,
    fixed_graphs,
    permutation_test,
    permute_trait_within_species,
    prepare_transfer,
    simulate_circular_boundary_world,
    split_species,
    transfer_statistic,
)


def test_null_permutation_preserves_coordinates_and_graph_geometry():
    world = simulate_circular_boundary_world(
        n_species=8,
        records_per_species=24,
        shared_fraction=0.0,
        amplitude=2.0,
        seed=4,
    )
    graphs = fixed_graphs(world.samples, k=3)
    rng = np.random.default_rng(99)
    for sample in world.samples:
        permuted = permute_trait_within_species(sample, rng)
        assert np.array_equal(permuted.coordinates, sample.coordinates)
        before = build_species_edges(sample, edge_nodes=graphs[sample.species])
        after = build_species_edges(permuted, edge_nodes=graphs[sample.species])
        assert np.array_equal(before.nodes, after.nodes)
        assert np.array_equal(before.midpoint, after.midpoint)


def test_prepared_geometry_matches_direct_edge_level_transfer():
    world = simulate_circular_boundary_world(
        n_species=10,
        records_per_species=28,
        shared_fraction=1.0,
        amplitude=2.0,
        noise_sd=0.5,
        transition_width=0.15,
        seed=11,
    )
    train, evaluation = split_species(
        [s.species for s in world.samples], eval_fraction=0.4, seed=3
    )
    edges = {s.species: build_species_edges(s, k=4) for s in world.samples}
    train_edges = [edges[s] for s in train]
    eval_edges = [edges[s] for s in evaluation]
    direct = transfer_statistic(
        train_edges,
        eval_edges,
        bandwidth=0.25,
        segment_points=5,
    )
    prepared = prepare_transfer(
        train_edges,
        eval_edges,
        bandwidth=0.25,
        segment_points=5,
    )
    fast = prepared.score(
        {s: edges[s].turnover for s in train},
        {s: edges[s].turnover for s in evaluation},
    )
    assert np.isclose(direct.statistic, fast.statistic, atol=1e-12)
    assert direct.species_scores.keys() == fast.species_scores.keys()


def test_full_shared_strong_world_beats_its_permutation_null():
    world = simulate_circular_boundary_world(
        n_species=16,
        records_per_species=40,
        shared_fraction=1.0,
        amplitude=3.0,
        noise_sd=0.5,
        transition_width=0.15,
        seed=101,
    )
    train, evaluation = split_species(
        [s.species for s in world.samples], eval_fraction=0.5, seed=17
    )
    result = permutation_test(
        world.samples,
        train_species=train,
        eval_species=evaluation,
        k=4,
        bandwidth=0.25,
        n_permutations=19,
        seed=202,
    )
    assert result.observed.statistic > result.null_mean + 0.05
    assert result.p_value <= 0.10


def test_zero_shared_strong_world_does_not_gain_from_other_species():
    world = simulate_circular_boundary_world(
        n_species=16,
        records_per_species=40,
        shared_fraction=0.0,
        amplitude=3.0,
        noise_sd=0.5,
        transition_width=0.15,
        seed=303,
    )
    train, evaluation = split_species(
        [s.species for s in world.samples], eval_fraction=0.5, seed=23
    )
    result = permutation_test(
        world.samples,
        train_species=train,
        eval_species=evaluation,
        k=4,
        bandwidth=0.25,
        n_permutations=19,
        seed=404,
    )
    assert abs(result.observed.statistic - result.null_mean) < 0.12
