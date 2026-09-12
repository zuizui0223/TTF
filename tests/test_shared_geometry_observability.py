import numpy as np

from ttf.core import knn_edges
from ttf.geometry import SpeciesGeometry
from ttf.shared_geometry_observability import shared_transition_observability


def _toy_geometry():
    a = SpeciesGeometry(
        species="a",
        coordinates=np.array(
            [[-2.0, -0.3], [-1.0, 0.2], [0.0, -0.1], [1.0, 0.1], [2.0, 0.0], [0.5, 0.7]],
            dtype=float,
        ),
    )
    b = SpeciesGeometry(
        species="b",
        coordinates=np.array(
            [[-1.5, 1.1], [-0.8, 0.8], [0.1, 1.0], [0.9, 0.9], [1.7, 1.2], [-0.2, 1.5]],
            dtype=float,
        ),
    )
    return (a, b)


def test_shared_observability_is_deterministic_and_bounded():
    geometries = _toy_geometry()
    graphs = {g.species: knn_edges(g.coordinates, k=2) for g in geometries}
    first = shared_transition_observability(
        geometries,
        graphs,
        transition_width=0.2,
        n_directions=64,
        seed=123,
    )
    second = shared_transition_observability(
        geometries,
        graphs,
        transition_width=0.2,
        n_directions=64,
        seed=123,
    )
    np.testing.assert_allclose(
        first.equal_species_mean_contrast_by_direction,
        second.equal_species_mean_contrast_by_direction,
        atol=0.0,
        rtol=0.0,
    )
    assert first.directions == 64
    assert np.all(first.equal_species_mean_contrast_by_direction >= 0.0)
    assert np.all(first.equal_species_mean_contrast_by_direction <= 2.0)
    assert set(first.species_mean_contrast) == {"a", "b"}
    assert all(0.0 <= value <= 2.0 for value in first.species_mean_contrast.values())
    assert all(0.0 <= value <= 2.0 for value in first.species_median_contrast.values())
    assert all(0.0 <= value <= 1.0 for value in first.species_nontrivial_fraction.values())


def test_shared_observability_changes_when_geometry_changes():
    geometries = list(_toy_geometry())
    graphs = {g.species: knn_edges(g.coordinates, k=2) for g in geometries}
    base = shared_transition_observability(
        geometries,
        graphs,
        n_directions=64,
        seed=321,
    )
    stretched = SpeciesGeometry(
        species="b",
        coordinates=geometries[1].coordinates * np.array([3.0, 1.0]),
    )
    changed_geometries = (geometries[0], stretched)
    changed_graphs = {
        g.species: knn_edges(g.coordinates, k=2) for g in changed_geometries
    }
    changed = shared_transition_observability(
        changed_geometries,
        changed_graphs,
        n_directions=64,
        seed=321,
    )
    assert not np.allclose(
        base.equal_species_mean_contrast_by_direction,
        changed.equal_species_mean_contrast_by_direction,
    )
