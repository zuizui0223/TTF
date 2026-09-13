import numpy as np

from ttf.genetic_geometry import prepare_density_scaled_genetic_geometry
from ttf.genetic_ibd import crossfit_ibd_residuals
from ttf.genetic_simulate import simulate_genetic_distance_world


def _geometry(n=24):
    theta = np.linspace(0.0, 2.0 * np.pi, n, endpoint=False)
    coordinates = np.column_stack([np.cos(theta), np.sin(theta)])
    return prepare_density_scaled_genetic_geometry(coordinates)


def _geographic(geometry):
    nodes = geometry.edge_nodes
    return np.linalg.norm(
        geometry.coordinates[nodes[:, 1]] - geometry.coordinates[nodes[:, 0]],
        axis=1,
    )


def test_ibd_only_noise_free_world_uses_exact_canonical_rank_scale():
    geometry = _geometry()
    world = simulate_genetic_distance_world(
        {"sp": geometry},
        shared_fraction=0.0,
        residual_amplitude=0.0,
        ibd_strength=2.5,
        noise_sd=0.0,
        seed=7,
    )
    geographic = _geographic(geometry)
    assert np.array_equal(world.genetic_distance["sp"], geographic)
    result = crossfit_ibd_residuals(
        world.genetic_distance["sp"],
        geographic,
        geometry.edge_nodes,
        min_training_edges=5,
    )
    assert np.array_equal(result.residual, np.zeros_like(result.residual))


def test_mixed_positive_ibd_strengths_share_same_exact_rank_canonicalization():
    geometries = {"a": _geometry(), "b": _geometry()}
    world = simulate_genetic_distance_world(
        geometries,
        shared_fraction=0.0,
        residual_amplitude=0.0,
        ibd_strength={"a": 0.5, "b": 3.0},
        noise_sd=0.0,
        seed=11,
    )
    assert np.array_equal(world.genetic_distance["a"], _geographic(geometries["a"]))
    assert np.array_equal(world.genetic_distance["b"], _geographic(geometries["b"]))
    assert np.max(world.residual_edge_signal["a"]) == 0.0
    assert np.max(world.residual_edge_signal["b"]) == 0.0


def test_zero_ibd_strength_remains_zero_in_pure_arm():
    geometry = _geometry()
    world = simulate_genetic_distance_world(
        {"sp": geometry},
        shared_fraction=0.0,
        residual_amplitude=0.0,
        ibd_strength=0.0,
        noise_sd=0.0,
        seed=12,
    )
    assert np.array_equal(
        world.genetic_distance["sp"],
        np.zeros(geometry.n_edges, dtype=float),
    )


def test_shared_fraction_controls_only_residual_field_membership():
    geometries = {f"sp{i}": _geometry() for i in range(8)}
    private = simulate_genetic_distance_world(
        geometries,
        shared_fraction=0.0,
        residual_amplitude=1.0,
        noise_sd=0.0,
        seed=13,
    )
    shared = simulate_genetic_distance_world(
        geometries,
        shared_fraction=1.0,
        residual_amplitude=1.0,
        noise_sd=0.0,
        seed=13,
    )
    assert private.shared_species == ()
    assert shared.shared_species == tuple(sorted(geometries))
