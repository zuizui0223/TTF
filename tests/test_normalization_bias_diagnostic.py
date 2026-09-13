import numpy as np

from ttf.conditional_null_centering import (
    identical_uniform_circle_samples,
    q51_private_relation_scores_for_phases,
)
from ttf.geometry_control import length_orthogonalized_turnover
from ttf.heterogeneous_simulate import simulate_q5_world
from ttf.normalization_bias_diagnostic import (
    length_residual_unscaled,
    normalization_bias_phase_redraw_scores,
    normalization_bias_scores_for_phases,
)


def _split():
    return (
        [f"sp_{i:03d}" for i in range(20)],
        [f"sp_{i:03d}" for i in range(20, 40)],
    )


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


def test_unscaled_residual_differs_from_q51_only_by_final_scale():
    rng = np.random.default_rng(771)
    turnover = rng.normal(size=41)
    length = np.exp(rng.normal(size=41))
    raw = length_residual_unscaled(turnover, length)
    normalized = length_orthogonalized_turnover(turnover, length)
    sd = float(np.std(raw, ddof=0))
    assert sd > np.sqrt(np.finfo(float).eps)
    n = len(raw)
    target_sd = float(np.sqrt((n * n - 1.0) / (12.0 * n * n)))
    assert np.allclose(
        normalized,
        raw * (target_sd / sd),
        atol=1e-12,
        rtol=0.0,
    )
    assert abs(float(np.mean(raw))) < 1e-14


def test_normalized_arm_exactly_reproduces_frozen_conditional_null_scorer():
    world = simulate_q5_world(
        response_mode="null_component",
        geometry_profile="matched",
        seed=93417,
    )
    names = _split()[0] + _split()[1]
    rng = np.random.default_rng(44119)
    phases = {name: rng.uniform(0.0, 2.0 * np.pi, size=5) for name in names}

    paired = normalization_bias_scores_for_phases(
        world.samples,
        phases_by_species=phases,
        **_kwargs(),
    )
    frozen = q51_private_relation_scores_for_phases(
        world.samples,
        phases_by_species=phases,
        **_kwargs(),
    )
    assert np.allclose(
        paired.normalized_statistics,
        frozen.statistics,
        atol=1e-12,
        rtol=0.0,
    )
    for name in _split()[1]:
        assert np.allclose(
            paired.normalized_species_scores[name],
            frozen.species_scores[name],
            atol=1e-12,
            rtol=0.0,
        )


def test_paired_arms_share_phases_graphs_and_are_reproducible():
    samples = identical_uniform_circle_samples(seed=818)
    first = normalization_bias_phase_redraw_scores(
        samples,
        n_phase_draws=7,
        phase_seed=20261,
        **_kwargs(),
    )
    second = normalization_bias_phase_redraw_scores(
        samples,
        n_phase_draws=7,
        phase_seed=20261,
        **_kwargs(),
    )
    assert first.normalized_statistics.shape == (7,)
    assert first.unscaled_statistics.shape == (7,)
    assert np.isfinite(first.normalized_statistics).all()
    assert np.isfinite(first.unscaled_statistics).all()
    assert np.allclose(first.normalized_statistics, second.normalized_statistics, atol=0.0, rtol=0.0)
    assert np.allclose(first.unscaled_statistics, second.unscaled_statistics, atol=0.0, rtol=0.0)
    assert first.graph_k == second.graph_k
    assert first.effective_n == second.effective_n
    for name in first.phases:
        assert np.array_equal(first.phases[name], second.phases[name])
