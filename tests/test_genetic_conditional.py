from __future__ import annotations

import numpy as np

from ttf.genetic_conditional import (
    prepare_genetic_conditional_design,
    score_genetic_conditional_world,
)
from ttf.genetic_conditional_simulate import (
    simulate_grouped_genetic_distance_world,
)
from ttf.genetic_geometry import (
    prepare_density_scaled_genetic_geometry,
)


def _geometries():
    theta = np.linspace(
        0,
        2 * np.pi,
        16,
        endpoint=False,
    )
    base = np.column_stack([
        np.cos(theta),
        np.sin(theta),
    ])
    out = {}
    for index in range(12):
        shift = np.array([
            0.2 * (index // 4),
            0.2 * (index // 4),
        ])
        out[f"sp{index:02d}"] = (
            prepare_density_scaled_genetic_geometry(
                base + shift
            )
        )
    return out


def test_grouped_world_and_conditional_design_are_outcome_independent() -> None:
    geometries = _geometries()
    names = tuple(sorted(geometries))
    train = names[:6]
    evaluation = names[6:]
    taxonomy = {
        name: {
            "order": (
                "A"
                if int(name[-2:]) % 2 == 0
                else "B"
            )
        }
        for name in names
    }
    design = prepare_genetic_conditional_design(
        geometries,
        taxonomy,
        train_species=train,
        eval_species=evaluation,
        bandwidth=2.0,
        radius=5.0,
        coverage_threshold=0.25,
        min_sources=2,
    )
    assert (
        design.source_sets.primary_eval_species
        == evaluation
    )
    assert (
        design.source_sets.secondary_eval_species
        == evaluation
    )

    groups = {
        name: taxonomy[name]["order"]
        for name in names
    }
    world = simulate_grouped_genetic_distance_world(
        geometries,
        groups,
        residual_amplitude=2.0,
        noise_sd=0.05,
        seed=42,
    )
    score = score_genetic_conditional_world(
        design,
        world.genetic_distance,
    )
    assert score.delta_geographic.shape == (1,)
    assert score.delta_order_given_geography.shape == (1,)
    assert np.isfinite(score.delta_geographic[0])
    assert np.isfinite(
        score.delta_order_given_geography[0]
    )


def test_private_group_world_uses_one_field_per_species() -> None:
    geometries = _geometries()
    groups = {
        name: name
        for name in geometries
    }
    world = simulate_grouped_genetic_distance_world(
        geometries,
        groups,
        residual_amplitude=1.0,
        noise_sd=0.1,
        seed=7,
    )
    assert len(world.group_normal) == len(geometries)
    assert set(world.group_by_species) == set(geometries)
