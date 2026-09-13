import numpy as np

from ttf.genetic_geometry import (
    endpoint_disjoint_training_counts,
    prepare_density_scaled_genetic_geometry,
    prepare_genetic_sampling_geometry,
)


def test_prepare_genetic_geometry_collapses_only_exact_duplicate_localities():
    coordinates = np.asarray(
        [
            [0.0, 0.0],
            [0.0, 0.0],
            [0.00001, 0.0],
            [1.0, 0.0],
            [2.0, 0.0],
            [3.0, 0.0],
            [4.0, 0.0],
            [5.0, 0.0],
        ]
    )
    geometry = prepare_genetic_sampling_geometry(coordinates, k=2)
    assert geometry.n_records == 8
    assert geometry.n_localities == 7
    assert geometry.graph_k == 2
    assert sorted(geometry.records_per_locality.tolist()) == [1, 1, 1, 1, 1, 1, 2]
    assert geometry.record_to_locality[0] == geometry.record_to_locality[1]
    assert geometry.record_to_locality[2] != geometry.record_to_locality[0]


def test_endpoint_disjoint_counts_depend_only_on_graph_endpoints():
    edges = np.asarray(
        [
            [0, 1],
            [0, 2],
            [1, 2],
            [2, 3],
            [3, 4],
            [4, 5],
        ],
        dtype=np.int64,
    )
    counts = endpoint_disjoint_training_counts(edges)
    assert counts.shape == (len(edges),)
    assert counts[0] == 3  # (2,3), (3,4), (4,5)
    assert counts[-1] == 4  # all except edges touching 4 or 5


def test_genetic_geometry_reports_crossfit_information_before_outcomes_exist():
    rng = np.random.default_rng(20260913)
    coordinates = rng.normal(size=(16, 2))
    geometry = prepare_genetic_sampling_geometry(coordinates, k=4)
    assert geometry.n_localities == 16
    assert geometry.n_edges > geometry.n_localities
    assert geometry.min_endpoint_disjoint_training_edges >= 3


def test_density_scaled_genetic_geometry_uses_unique_locality_count_only():
    x = np.arange(20.0)
    coordinates = np.column_stack([x, np.sin(x)])
    coordinates = np.vstack([coordinates, coordinates[:5]])
    geometry = prepare_density_scaled_genetic_geometry(
        coordinates,
        neighbor_fraction=0.15,
    )
    assert geometry.n_records == 25
    assert geometry.n_localities == 20
    assert geometry.graph_k == 3  # floor(0.15 * 20 + 0.5)
