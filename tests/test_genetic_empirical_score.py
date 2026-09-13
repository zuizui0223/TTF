from __future__ import annotations

import numpy as np

from ttf.genetic_empirical_score import (
    score_genetic_distance_mapping,
    score_self_detectability_mapping,
)
from ttf.genetic_gate import prepare_genetic_ttf_design, score_genetic_world
from ttf.genetic_geometry import prepare_density_scaled_genetic_geometry
from ttf.genetic_self_detectability import (
    prepare_genetic_self_detectability,
    score_genetic_self_world_batch,
)
from ttf.genetic_simulate import simulate_genetic_distance_world


def _geometries() -> dict:
    out = {}
    for species_index in range(8):
        coordinates = np.column_stack(
            [
                np.arange(12, dtype=float) * 25.0 + species_index * 3.0,
                np.sin(np.arange(12, dtype=float) / 2.0) * 10.0,
                np.full(12, species_index * 2.0, dtype=float),
            ]
        )
        out[f"species_{species_index}"] = prepare_density_scaled_genetic_geometry(
            coordinates, neighbor_fraction=0.15
        )
    return out


def test_empirical_cross_adapter_matches_synthetic_scorer_exactly() -> None:
    geometries = _geometries()
    names = tuple(sorted(geometries))
    design = prepare_genetic_ttf_design(
        geometries,
        train_species=names[:4],
        eval_species=names[4:],
        bandwidth=500.0,
        prior_strength=0.25,
        segment_points=5,
        min_training_edges=5,
        strength_neighbours=4,
    )
    world = simulate_genetic_distance_world(
        geometries,
        shared_fraction=0.5,
        residual_amplitude=2.0,
        ibd_strength=1.0,
        noise_sd=0.1,
        transition_width=0.2,
        noise_dimensions=2,
        seed=12345,
    )
    expected = score_genetic_world(design, world)
    observed = score_genetic_distance_mapping(design, world.genetic_distance)
    assert observed.statistic == expected.statistic
    assert observed.training_strength == expected.training_strength
    assert observed.species_scores == expected.species_scores
    for name in names:
        np.testing.assert_array_equal(
            observed.ibd[name].residual_turnover,
            expected.ibd[name].residual_turnover,
        )


def test_empirical_self_adapter_matches_synthetic_self_scorer_exactly() -> None:
    geometries = _geometries()
    names = tuple(sorted(geometries))
    design = prepare_genetic_ttf_design(
        geometries,
        train_species=names[:4],
        eval_species=names[4:],
        bandwidth=500.0,
        prior_strength=0.25,
        segment_points=5,
        min_training_edges=5,
        strength_neighbours=4,
    )
    self_design = prepare_genetic_self_detectability(
        design,
        bandwidth=500.0,
        prior_strength=0.25,
        prior_mean=0.5,
        segment_points=5,
    )
    world = simulate_genetic_distance_world(
        geometries,
        shared_fraction=0.0,
        residual_amplitude=2.0,
        ibd_strength=1.0,
        noise_sd=0.1,
        transition_width=0.2,
        noise_dimensions=2,
        seed=67890,
    )
    expected = score_genetic_self_world_batch(design, self_design, [world])
    observed = score_self_detectability_mapping(
        design, self_design, world.genetic_distance
    )
    assert observed.statistic == float(expected.statistics[0])
    assert observed.species_scores == {
        name: float(expected.species_scores[name][0]) for name in self_design.species
    }
