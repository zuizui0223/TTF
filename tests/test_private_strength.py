import numpy as np

from ttf.private_strength import (
    edge_midpoint_neighbor_coherence,
    edge_midpoint_neighbor_indices,
    training_private_strength,
    training_private_strength_from_indices,
)


def test_neighbor_coherence_is_high_for_local_structure_and_low_for_alternation():
    midpoint = np.arange(12, dtype=float)[:, None]
    smooth = np.r_[np.zeros(6), np.ones(6)]
    alternating = np.arange(12, dtype=float) % 2
    assert edge_midpoint_neighbor_coherence(smooth, midpoint, k=2) > 0.7
    assert edge_midpoint_neighbor_coherence(alternating, midpoint, k=2) < 0.0


def test_training_strength_is_equal_species_mean_and_deterministic():
    midpoint = {
        "a": np.arange(8, dtype=float)[:, None],
        "b": np.arange(8, dtype=float)[:, None],
    }
    responses = {
        "a": np.r_[np.zeros(4), np.ones(4)],
        "b": np.linspace(-1.0, 1.0, 8),
    }
    value, per_species = training_private_strength(
        responses,
        midpoint,
        ["a", "b"],
        k=2,
    )
    expected = np.mean(list(per_species.values()))
    assert value == expected
    again, again_species = training_private_strength(
        responses,
        midpoint,
        ["a", "b"],
        k=2,
    )
    assert value == again
    assert per_species == again_species


def test_precomputed_neighbour_fast_path_is_exact():
    midpoint = {
        "a": np.column_stack([np.arange(9, dtype=float), np.zeros(9)]),
        "b": np.column_stack([np.arange(9, dtype=float), np.ones(9)]),
    }
    responses = {
        "a": np.linspace(-0.4, 0.6, 9),
        "b": np.r_[np.zeros(4), np.ones(5)],
    }
    direct = training_private_strength(
        responses, midpoint, ["a", "b"], k=3
    )
    indices = {
        name: edge_midpoint_neighbor_indices(points, k=3)
        for name, points in midpoint.items()
    }
    fast = training_private_strength_from_indices(
        responses, indices, ["a", "b"]
    )
    assert direct[0] == fast[0]
    assert direct[1] == fast[1]
