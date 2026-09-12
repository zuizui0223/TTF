import numpy as np

from ttf.conditional_null_centering import (
    identical_uniform_circle_samples,
    q51_private_relation_phase_redraw_scores,
    q51_private_relation_scores_for_phases,
)
from ttf.heterogeneous_simulate import simulate_q5_world
from ttf.mismatch import PairedSpeciesSample
from ttf.q51_inference import q51_paired_heldout_species_bootstrap_test


def _split():
    train = [f"sp_{i:03d}" for i in range(20)]
    evaluation = [f"sp_{i:03d}" for i in range(20, 40)]
    return train, evaluation


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


def test_one_phase_draw_reproduces_frozen_q51_coupling_statistic():
    world = simulate_q5_world(
        response_mode="null_private_relation",
        geometry_profile="matched",
        amplitude=2.0,
        noise_sd=0.8,
        transition_width=0.2,
        seed=91221,
    )
    phases = {
        name: np.asarray([float(world.private_phase[name])], dtype=float)
        for name in world.private_phase
    }
    direct = q51_private_relation_scores_for_phases(
        world.samples,
        phases_by_species=phases,
        **_kwargs(),
    )
    frozen = q51_paired_heldout_species_bootstrap_test(
        world.samples,
        **_kwargs(),
        n_bootstrap=19,
        seed=991,
    )
    assert direct.statistics.shape == (1,)
    assert np.isclose(
        direct.statistics[0],
        frozen.coupling_statistic,
        atol=1e-12,
        rtol=0.0,
    )
    for name in _split()[1]:
        assert np.isclose(
            direct.species_scores[name][0],
            frozen.coupling_species_scores[name],
            atol=1e-12,
            rtol=0.0,
        )


def test_phase_redraw_scores_ignore_input_response_values_on_fixed_geometry():
    samples = identical_uniform_circle_samples(seed=73)
    changed = tuple(
        PairedSpeciesSample(
            species=sample.species,
            coordinates=sample.coordinates,
            state_a=np.linspace(-3.0, 4.0, len(sample.coordinates)) + i,
            state_b=np.cos(np.arange(len(sample.coordinates), dtype=float) + i),
        )
        for i, sample in enumerate(samples)
    )
    first = q51_private_relation_phase_redraw_scores(
        samples,
        n_phase_draws=8,
        phase_seed=5521,
        **_kwargs(),
    )
    second = q51_private_relation_phase_redraw_scores(
        changed,
        n_phase_draws=8,
        phase_seed=5521,
        **_kwargs(),
    )
    assert np.allclose(first.statistics, second.statistics, atol=1e-12, rtol=0.0)
    assert first.graph_k == second.graph_k
    assert first.effective_n == second.effective_n


def test_phase_redraw_scores_are_finite_and_reproducible():
    world = simulate_q5_world(
        response_mode="null_component",
        geometry_profile="shifted",
        seed=33118,
    )
    first = q51_private_relation_phase_redraw_scores(
        world.samples,
        n_phase_draws=12,
        phase_seed=8803,
        **_kwargs(),
    )
    second = q51_private_relation_phase_redraw_scores(
        world.samples,
        n_phase_draws=12,
        phase_seed=8803,
        **_kwargs(),
    )
    assert first.statistics.shape == (12,)
    assert np.isfinite(first.statistics).all()
    assert np.allclose(first.statistics, second.statistics, atol=0.0, rtol=0.0)
