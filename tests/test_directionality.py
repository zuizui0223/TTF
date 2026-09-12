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


def test_strict_private_periodic_world_separates_absolute_space_from_local_direction():
    world = simulate_directional_world(
        mode="strict_private_periodic_directional_change",
        n_species=20,
        records_per_species=60,
        noise_sd=0.0,
        strict_private_arc_halfwidth=1.0,
        seed=41,
    )
    phases = np.asarray(list(world.front_by_species.values()), dtype=float)
    assert np.all((phases >= 0.0) & (phases < 2.0 * np.pi))
    assert np.ptp(phases) > np.pi
    assert all(x.sample.coordinates.shape == (60, 2) for x in world.samples)
    assert all(np.min(x.orientation) < 0.0 < np.max(x.orientation) for x in world.samples)

    # Directional mismatch is common in each system's predeclared local
    # orientation even though the absolute spatial phase differs by system.
    labels = [x.species for x in world.samples]
    deltas = heldout_oriented_mismatch_deltas(
        world.samples,
        eval_species=labels,
        front_center=0.0,
        front_window=0.35,
    )
    vec = np.asarray([deltas[s] for s in labels if s in deltas], dtype=float)
    assert len(vec) >= 18
    assert np.mean(vec > 0.0) >= 0.9
    assert vec.mean() > 1.0


def test_shared_front_zone_is_old_bounded_private_world_with_new_semantics_only():
    old = simulate_directional_world(
        mode="private_directional_change",
        n_species=12,
        records_per_species=40,
        noise_sd=0.25,
        seed=73,
    )
    zone = simulate_directional_world(
        mode="shared_front_zone_directional_change",
        n_species=12,
        records_per_species=40,
        noise_sd=0.25,
        seed=73,
    )
    assert old.front_by_species == zone.front_by_species
    for a, b in zip(old.samples, zone.samples, strict=True):
        np.testing.assert_array_equal(a.sample.coordinates, b.sample.coordinates)
        np.testing.assert_array_equal(a.orientation, b.orientation)
        np.testing.assert_array_equal(a.sample.state_a, b.sample.state_a)
        np.testing.assert_array_equal(a.sample.state_b, b.sample.state_b)
