import numpy as np

from ttf.private_strength import (
    edge_midpoint_neighbor_coherence,
    training_private_strength,
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
