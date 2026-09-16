from __future__ import annotations

import numpy as np
import pytest

from ttf.core import knn_edges
from ttf.genetic_ibd import (
    crossfit_ibd_residuals_prepared_reference,
    prepare_crossfit_ibd_design,
)
from ttf.phylogatr_compact_ibd import (
    crossfit_ibd_residuals_compact,
    prepare_compact_crossfit_ibd_design,
)


def _case(*, ties: bool) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    rng = np.random.default_rng(20260916)
    coordinates = rng.normal(size=(18, 3))
    nodes = knn_edges(coordinates, k=4)
    geographic = np.linalg.norm(
        coordinates[nodes[:, 0]] - coordinates[nodes[:, 1]], axis=1
    )
    genetic = np.square(rng.normal(size=len(nodes)))
    if ties:
        geographic = np.round(geographic, 1)
        genetic = np.round(genetic, 1)
    return geographic, nodes, genetic


@pytest.mark.parametrize("ties", [False, True])
def test_compact_ibd_matches_reference_with_and_without_ties(ties: bool) -> None:
    geographic, nodes, genetic = _case(ties=ties)
    reference_design = prepare_crossfit_ibd_design(
        geographic,
        nodes,
        min_training_edges=5,
    )
    reference = crossfit_ibd_residuals_prepared_reference(
        genetic,
        reference_design,
    )

    compact_design = prepare_compact_crossfit_ibd_design(
        geographic,
        nodes,
        min_training_edges=5,
    )
    compact = crossfit_ibd_residuals_compact(genetic, compact_design)

    assert np.allclose(compact.residual, reference.residual, rtol=0.0, atol=5e-14)
    assert np.array_equal(compact.residual_turnover, reference.residual_turnover)
    assert np.allclose(
        compact.observed_rank_fraction,
        reference.observed_rank_fraction,
        rtol=0.0,
        atol=5e-14,
    )
    assert np.allclose(
        compact.expected_rank_fraction,
        reference.expected_rank_fraction,
        rtol=0.0,
        atol=5e-14,
    )
    assert np.allclose(
        compact.geographic_rank_fraction,
        reference.geographic_rank_fraction,
        rtol=0.0,
        atol=5e-14,
    )
    assert np.array_equal(compact.n_training_edges, reference.n_training_edges)


def test_compact_design_does_not_store_edge_by_edge_training_arrays() -> None:
    geographic, nodes, _ = _case(ties=False)
    compact = prepare_compact_crossfit_ibd_design(
        geographic,
        nodes,
        min_training_edges=5,
    )

    assert not hasattr(compact, "training_indices")
    assert not hasattr(compact, "geographic_centered_ranks")
    assert compact.vertex_rank_correction.shape == (
        int(np.max(nodes)) + 1,
        len(nodes),
    )
    assert compact.vertex_rank_correction.nbytes < len(nodes) * len(nodes) * 8
