import numpy as np

from ttf.directionality import (
    directional_coupling_test,
    estimate_training_front_center,
    heldout_oriented_mismatch_deltas,
)
from ttf.directionality_simulate import simulate_directional_world


def _split(world):
    labels = [x.species for x in world.samples]
    half = len(labels) // 2
    return labels[:half], labels[half:]


def _deltas(mode: str, *, noise_sd: float = 0.0):
    world = simulate_directional_world(
        mode=mode,
        n_species=20,
        records_per_species=60,
        noise_sd=noise_sd,
        seed=13,
    )
    train, evaluation = _split(world)
    center = estimate_training_front_center(
        world.samples,
        train_species=train,
        k=9,
        front_edge_fraction=0.15,
    )
    deltas = heldout_oriented_mismatch_deltas(
        world.samples,
        eval_species=evaluation,
        front_center=center,
        front_window=0.35,
    )
    return center, np.asarray([deltas[s] for s in evaluation if s in deltas], dtype=float)


def test_shared_breakdown_and_recoupling_have_opposite_oriented_signs():
    c_break, d_break = _deltas("shared_breakdown")
    c_recover, d_recover = _deltas("shared_recoupling")
    assert abs(c_break) < 0.2
    assert abs(c_recover) < 0.2
    assert d_break.mean() > 1.0
    assert np.mean(d_break > 0) >= 0.9
    assert d_recover.mean() < -1.0
    assert np.mean(d_recover < 0) >= 0.9


def test_neutral_relation_rotation_has_zero_mismatch_direction():
    center, delta = _deltas("shared_neutral_rotation")
    assert abs(center) < 0.2
    assert np.max(np.abs(delta)) < 1e-12


def test_sign_flip_cancels_shared_direction_across_systems():
    center, delta = _deltas("shared_sign_flip")
    assert abs(center) < 0.2
    assert abs(float(delta.mean())) < 0.4
    assert np.mean(delta > 0) >= 0.3
    assert np.mean(delta < 0) >= 0.3


def test_directional_coupling_requires_front_and_heldout_sign():
    world = simulate_directional_world(
        mode="shared_breakdown",
        n_species=20,
        records_per_species=60,
        noise_sd=0.0,
        seed=29,
    )
    train, evaluation = _split(world)
    result = directional_coupling_test(
        world.samples,
        train_species=train,
        eval_species=evaluation,
        k=9,
        bandwidth=0.2,
        front_edge_fraction=0.15,
        front_window=0.35,
        sign_consistency_fraction=0.75,
        alpha=0.05,
        n_bootstrap=99,
        seed=31,
    )
    assert result.relation_detected
    assert result.breakdown_label
    assert not result.recoupling_label
    assert result.positive_fraction >= 0.9
