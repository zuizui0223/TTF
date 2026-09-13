import numpy as np

from ttf.genetic_geometry import prepare_density_scaled_genetic_geometry
from ttf.genetic_ibd import crossfit_ibd_residuals
from ttf.genetic_simulate import simulate_genetic_distance_world


def _geometry(n=24):
    theta = np.linspace(0.0, 2.0 * np.pi, n, endpoint=False)
    coordinates = np.column_stack([np.cos(theta), np.sin(theta)])
    return prepare_density_scaled_genetic_geometry(coordinates)


def test_ibd_only_noise_free_world_is_exactly_removed_by_crossfit():
    geometry = _geometry()
    world = simulate_genetic_distance_world(
        {"sp": geometry},
        shared_fraction=0.0,
        residual_amplitude=0.0,
        ibd_strength=2.5,
        noise_sd=0.0,
        seed=7,
    )
    nodes = geometry.edge_nodes
    start = geometry.coordinates[nodes[:, 0]]
    end = geometry.coordinates[nodes[:, 1]]
    geographic = np.linalg.norm(end - start, axis=1)
    result = crossfit_ibd_residuals(
        world.genetic_distance["sp"],
        geographic,
        nodes,
        min_training_edges=5,
    )
    assert np.allclose(result.residual, 0.0, atol=1e-12, rtol=0.0)


def test_mixed_ibd_strengths_are_supported_without_changing_geometry():
    geometries = {"a": _geometry(), "b": _geometry()}
    world = simulate_genetic_distance_world(
        geometries,
        shared_fraction=0.0,
        residual_amplitude=0.0,
        ibd_strength={"a": 0.5, "b": 3.0},
        noise_sd=0.0,
        seed=11,
    )
    assert set(world.genetic_distance) == {"a", "b"}
    assert len(world.genetic_distance["a"]) == geometries["a"].n_edges
    assert len(world.genetic_distance["b"]) == geometries["b"].n_edges
    assert np.max(world.residual_edge_signal["a"]) == 0.0
    assert np.max(world.residual_edge_signal["b"]) == 0.0


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
