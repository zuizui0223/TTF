from __future__ import annotations

import numpy as np

from ttf.genetic_conditional_simulate import simulate_grouped_genetic_distance_world
from ttf.genetic_geometry import prepare_density_scaled_genetic_geometry


def _geometries() -> dict[str, object]:
    theta = np.linspace(0.0, 2.0 * np.pi, 16, endpoint=False)
    base = np.column_stack([np.cos(theta), np.sin(theta)])
    return {
        f"sp{i:02d}": prepare_density_scaled_genetic_geometry(
            base + np.array([0.15 * (i // 4), 0.15 * (i // 4)])
        )
        for i in range(12)
    }


def test_grouped_world_is_reproducible_and_shares_one_field_per_group() -> None:
    geometries = _geometries()
    names = tuple(sorted(geometries))
    groups = {name: ("A" if int(name[-2:]) % 2 == 0 else "B") for name in names}
    first = simulate_grouped_genetic_distance_world(
        geometries,
        groups,
        residual_amplitude=2.0,
        noise_sd=0.1,
        seed=42,
    )
    second = simulate_grouped_genetic_distance_world(
        geometries,
        groups,
        residual_amplitude=2.0,
        noise_sd=0.1,
        seed=42,
    )
    assert set(first.group_normal) == {"A", "B"}
    assert set(first.group_offset) == {"A", "B"}
    for name in names:
        assert np.array_equal(first.genetic_distance[name], second.genetic_distance[name])
        assert np.array_equal(first.residual_edge_signal[name], second.residual_edge_signal[name])
    assert np.array_equal(first.group_normal["A"], second.group_normal["A"])


def test_private_null_is_unique_group_per_species() -> None:
    geometries = _geometries()
    groups = {name: name for name in geometries}
    world = simulate_grouped_genetic_distance_world(
        geometries,
        groups,
        residual_amplitude=2.0,
        noise_sd=0.1,
        seed=7,
    )
    assert len(world.group_normal) == len(geometries)
    assert len(world.group_offset) == len(geometries)
    assert world.group_by_species == groups


def test_group_labels_are_required_for_every_species_and_nonempty() -> None:
    geometries = _geometries()
    names = tuple(sorted(geometries))
    bad = {name: "A" for name in names[:-1]}
    try:
        simulate_grouped_genetic_distance_world(
            geometries, bad, residual_amplitude=1.0
        )
    except ValueError as exc:
        assert "every species" in str(exc)
    else:
        raise AssertionError("missing group label was accepted")

    blank = {name: "A" for name in names}
    blank[names[0]] = ""
    try:
        simulate_grouped_genetic_distance_world(
            geometries, blank, residual_amplitude=1.0
        )
    except ValueError as exc:
        assert "non-empty" in str(exc)
    else:
        raise AssertionError("blank group label was accepted")
