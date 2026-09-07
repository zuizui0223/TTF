import numpy as np
import pytest

from ttf.genetics import (
    PairwiseDistanceSample,
    build_pairwise_distance_edges,
    pairwise_sequence_distance_matrix,
    sequence_p_distance,
)


def test_sequence_p_distance_ignores_ambiguous_sites():
    assert np.isclose(sequence_p_distance("ACGTN-", "ACGAA?"), 1.0 / 4.0)


def test_pairwise_sequence_distance_matrix_is_symmetric_metric_input():
    matrix = pairwise_sequence_distance_matrix(["AAAA", "AAAT", "AATT"])
    assert np.allclose(matrix, matrix.T)
    assert np.allclose(np.diag(matrix), 0.0)
    assert np.isclose(matrix[0, 1], 0.25)
    assert np.isclose(matrix[0, 2], 0.50)


def test_pairwise_distance_sample_rejects_asymmetric_or_nonzero_diagonal():
    coords = np.array([[0.0, 0.0], [1.0, 0.0]])
    with pytest.raises(ValueError, match="symmetric"):
        PairwiseDistanceSample(
            "sp",
            coords,
            np.array([[0.0, 0.1], [0.2, 0.0]]),
        )
    with pytest.raises(ValueError, match="diagonal"):
        PairwiseDistanceSample(
            "sp",
            coords,
            np.array([[0.1, 0.2], [0.2, 0.0]]),
        )


def test_genetic_edges_read_fixed_matrix_without_modifying_it():
    coords = np.array(
        [[0.0, 0.0], [1.0, 0.0], [2.0, 0.0], [3.0, 0.0]], dtype=float
    )
    matrix = np.array(
        [
            [0.0, 0.1, 0.7, 0.8],
            [0.1, 0.0, 0.2, 0.9],
            [0.7, 0.2, 0.0, 0.6],
            [0.8, 0.9, 0.6, 0.0],
        ]
    )
    frozen = matrix.copy()
    sample = PairwiseDistanceSample("sp", coords, matrix)
    nodes = np.array([[0, 1], [1, 2], [2, 3]])
    edges = build_pairwise_distance_edges(sample, edge_nodes=nodes)
    assert np.array_equal(sample.distances, frozen)
    assert np.array_equal(edges.nodes, nodes)
    assert np.allclose(edges.turnover, [1.0 / 6.0, 3.0 / 6.0, 5.0 / 6.0])
