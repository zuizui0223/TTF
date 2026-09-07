import numpy as np

from ttf.semisynthetic import geometry_bandwidth, normalize_geometry, simulate_on_geometry


def toy_geometry():
    return {
        f"sp{i}": np.column_stack(
            [
                np.linspace(-1.0, 1.0, 8) + 0.05 * i,
                np.sin(np.linspace(-1.0, 1.0, 8) * np.pi) + 0.1 * i,
            ]
        )
        for i in range(8)
    }


def test_geometry_normalization_has_unit_pooled_rms_radius():
    normalized = normalize_geometry(toy_geometry())
    pooled = np.vstack(list(normalized.coordinates.values()))
    assert np.allclose(pooled.mean(axis=0), 0.0, atol=1e-12)
    rms = np.sqrt(np.mean(np.sum(pooled * pooled, axis=1)))
    assert np.isclose(rms, 1.0, atol=1e-12)


def test_geometry_bandwidth_is_positive_and_deterministic():
    g = toy_geometry()
    assert geometry_bandwidth(g, k=3) > 0
    assert geometry_bandwidth(g, k=3) == geometry_bandwidth(g, k=3)


def test_zero_shared_world_has_no_shared_species():
    world = simulate_on_geometry(
        toy_geometry(),
        shared_fraction=0.0,
        amplitude=2.0,
        min_shared_crossing_species=6,
        seed=11,
    )
    assert world.shared_species == ()
    assert world.shared_normal is None
    assert world.shared_offset is None
    assert len(world.samples) == 8


def test_full_shared_world_uses_one_boundary_crossed_by_minimum_species():
    world = simulate_on_geometry(
        toy_geometry(),
        shared_fraction=1.0,
        amplitude=2.0,
        min_shared_crossing_species=6,
        seed=19,
    )
    assert len(world.shared_species) == 8
    assert world.shared_normal is not None
    assert world.shared_offset is not None
    assert len(world.crossing_species) >= 6
    assert world.bandwidth > 0
