import numpy as np

from ttf.conditional_null_centering import _phase_turnover_matrix
from ttf.core import SpeciesEdges
from ttf.heterogeneous_simulate import simulate_q5_world
from ttf.mismatch import PairedSpeciesSample
from ttf.normalization_bias_diagnostic import (
    length_residual_unscaled,
    normalization_bias_scores_for_phases,
)
from ttf.phase_propensity_diagnostic import (
    edge_private_phase_crossing_probability,
    phase_expected_centered_rank_turnover,
    phase_expected_unscaled_q51_residual,
    phase_propensity_scores_for_phases,
)


def _manual_sample_and_edges():
    theta = np.asarray([0.0, 0.2, 0.7, 1.8, 3.0], dtype=float)
    coords = np.column_stack((np.cos(theta), np.sin(theta)))
    sample = PairedSpeciesSample(
        species="manual",
        coordinates=coords,
        state_a=np.zeros(len(theta)),
        state_b=np.zeros(len(theta)),
    )
    nodes = np.asarray([[0, 1], [0, 2], [0, 3], [1, 3], [2, 4], [0, 4]], dtype=np.int64)
    start = coords[nodes[:, 0]]
    end = coords[nodes[:, 1]]
    edges = SpeciesEdges(
        species="manual",
        nodes=nodes,
        start=start,
        end=end,
        midpoint=0.5 * (start + end),
        length=np.linalg.norm(end - start, axis=1),
        turnover=np.full(len(nodes), 0.5),
    )
    return sample, edges


def _split():
    return [f"sp_{i:03d}" for i in range(20)], [f"sp_{i:03d}" for i in range(20, 40)]


def _kwargs():
    train, evaluation = _split()
    return dict(
        train_species=train,
        eval_species=evaluation,
        graph_fraction=0.15,
        bandwidth=0.2,
        prior_strength=0.25,
        segment_points=5,
        edge_chunk_size=32,
        train_chunk_size=2048,
    )


def test_crossing_probability_matches_dense_uniform_phase_frequency():
    sample, edges = _manual_sample_and_edges()
    analytic = edge_private_phase_crossing_probability(sample, edges)
    phases = (np.arange(32768, dtype=float) + 0.37) * (2.0 * np.pi / 32768.0)
    theta = np.mod(np.arctan2(sample.coordinates[:, 1], sample.coordinates[:, 0]), 2.0 * np.pi)
    side = np.sin(theta[:, None] - phases[None, :]) >= 0.0
    observed = np.mean(
        side[edges.nodes[:, 0], :] != side[edges.nodes[:, 1], :],
        axis=1,
    )
    assert np.allclose(observed, analytic, atol=8e-5, rtol=0.0)


def test_analytic_phase_mean_centered_rank_matches_dense_phase_average():
    sample, edges = _manual_sample_and_edges()
    phases = (np.arange(16384, dtype=float) + 0.23) * (2.0 * np.pi / 16384.0)
    turnover = _phase_turnover_matrix(sample, edges, phases)
    observed = np.mean(turnover - 0.5, axis=1)
    analytic = phase_expected_centered_rank_turnover(sample, edges)
    assert np.allclose(observed, analytic, atol=2e-4, rtol=0.0)
    assert abs(float(np.mean(analytic))) < 1e-14


def test_analytic_phase_mean_unscaled_residual_matches_dense_phase_average():
    sample, edges = _manual_sample_and_edges()
    phases = (np.arange(16384, dtype=float) + 0.41) * (2.0 * np.pi / 16384.0)
    turnover = _phase_turnover_matrix(sample, edges, phases)
    observed = np.mean(
        np.column_stack(
            [length_residual_unscaled(turnover[:, j], edges.length) for j in range(turnover.shape[1])]
        ),
        axis=1,
    )
    analytic = phase_expected_unscaled_q51_residual(sample, edges)
    assert np.allclose(observed, analytic, atol=2e-4, rtol=0.0)


def test_phase_propensity_baseline_exactly_reproduces_frozen_unscaled_arm():
    world = simulate_q5_world(
        response_mode="null_component",
        geometry_profile="matched",
        seed=92431,
    )
    train, evaluation = _split()
    rng = np.random.default_rng(924311)
    phases = {
        name: rng.uniform(0.0, 2.0 * np.pi, size=3)
        for name in train + evaluation
    }
    existing = normalization_bias_scores_for_phases(
        world.samples,
        phases_by_species=phases,
        **_kwargs(),
    )
    diagnostic = phase_propensity_scores_for_phases(
        world.samples,
        phases_by_species=phases,
        **_kwargs(),
    )
    assert np.allclose(
        diagnostic.baseline_statistics,
        existing.unscaled_statistics,
        atol=1e-12,
        rtol=0.0,
    )
    for name in evaluation:
        assert np.allclose(
            diagnostic.baseline_species_scores[name],
            existing.unscaled_species_scores[name],
            atol=1e-12,
            rtol=0.0,
        )


def test_expected_training_residual_is_geometry_only():
    world = simulate_q5_world(
        response_mode="null_component",
        geometry_profile="shifted",
        seed=92432,
    )
    changed = tuple(
        PairedSpeciesSample(
            species=s.species,
            coordinates=s.coordinates,
            state_a=np.arange(len(s.coordinates), dtype=float),
            state_b=-np.arange(len(s.coordinates), dtype=float),
        )
        for s in world.samples
    )
    train, evaluation = _split()
    rng = np.random.default_rng(924322)
    phases = {
        name: rng.uniform(0.0, 2.0 * np.pi, size=2)
        for name in train + evaluation
    }
    first = phase_propensity_scores_for_phases(world.samples, phases_by_species=phases, **_kwargs())
    second = phase_propensity_scores_for_phases(changed, phases_by_species=phases, **_kwargs())
    for name in train:
        assert np.allclose(first.expected_training_residual[name], second.expected_training_residual[name], atol=0.0, rtol=0.0)
