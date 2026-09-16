from __future__ import annotations

import numpy as np
import pytest

from ttf.core import knn_edges
from ttf.genetic_batch_execution import (
    prepare_genetic_cached_transfer,
    score_genetic_world_batch,
)
from ttf.genetic_gate import prepare_genetic_ttf_design
from ttf.genetic_geometry import prepare_density_scaled_genetic_geometry
from ttf.genetic_ibd import (
    crossfit_ibd_residuals_prepared_reference,
    prepare_crossfit_ibd_design,
)
from ttf.genetic_simulate import simulate_genetic_distance_world
from ttf.phylogatr_compact_execution import (
    prepare_phylogatr_compact_cached_transfer,
    prepare_phylogatr_compact_ttf_design,
    score_phylogatr_compact_world_batch,
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


def _small_geometries() -> dict[str, object]:
    rng = np.random.default_rng(4147)
    out = {}
    for index in range(12):
        coordinates = rng.normal(size=(18, 3)) + np.array([index * 0.15, 0.0, 0.0])
        out[f"species_{index:02d}"] = prepare_density_scaled_genetic_geometry(
            coordinates,
            neighbor_fraction=0.15,
        )
    return out


def test_compact_phase3_batch_matches_existing_batch_end_to_end() -> None:
    geometries = _small_geometries()
    names = tuple(sorted(geometries))
    train = names[:6]
    evaluation = names[6:]
    kwargs = dict(
        train_species=train,
        eval_species=evaluation,
        bandwidth=500.0,
        prior_strength=0.25,
        segment_points=5,
        min_training_edges=5,
        strength_neighbours=4,
    )
    reference_design = prepare_genetic_ttf_design(geometries, **kwargs)
    compact_design = prepare_phylogatr_compact_ttf_design(geometries, **kwargs)

    worlds = tuple(
        simulate_genetic_distance_world(
            geometries,
            shared_fraction=shared,
            residual_amplitude=amplitude,
            ibd_strength=1.0,
            noise_sd=0.10,
            transition_width=0.20,
            noise_dimensions=2,
            seed=7000 + index,
        )
        for index, (shared, amplitude) in enumerate(((0.0, 1.0), (0.5, 2.0), (1.0, 2.0)))
    )

    reference_cached = prepare_genetic_cached_transfer(reference_design)
    compact_cached = prepare_phylogatr_compact_cached_transfer(compact_design)
    reference = score_genetic_world_batch(reference_design, worlds, reference_cached)
    compact = score_phylogatr_compact_world_batch(compact_design, worlds, compact_cached)

    assert np.allclose(compact.statistics, reference.statistics, rtol=0.0, atol=2e-13)
    assert np.allclose(
        compact.training_strengths,
        reference.training_strengths,
        rtol=0.0,
        atol=2e-13,
    )
    assert np.array_equal(compact.n_eval_species, reference.n_eval_species)
    for name in evaluation:
        assert np.allclose(
            compact.species_scores[name],
            reference.species_scores[name],
            rtol=0.0,
            atol=2e-13,
        )
