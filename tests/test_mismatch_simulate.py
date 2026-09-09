import numpy as np

from ttf.core import split_species
from ttf.mismatch import build_coupling_edges, build_mismatch_edges, pointwise_mismatch
from ttf.mismatch_simulate import simulate_paired_transition_world
from ttf.transfer import transfer_statistic


def test_q0_coupled_transition_has_no_mismatch_or_coupling_signal() -> None:
    world = simulate_paired_transition_world(
        mode="coupled_shared_transition",
        n_species=6,
        records_per_species=30,
        amplitude=2.0,
        noise_sd=0.2,
        seed=1,
    )
    for sample in world.samples:
        assert np.allclose(pointwise_mismatch(sample), 0.0)
        mismatch_edges = build_mismatch_edges(sample, k=4)
        coupling_edges = build_coupling_edges(sample, k=4)
        assert np.allclose(mismatch_edges.turnover, 0.5)
        assert np.allclose(coupling_edges.turnover, 0.5)


def test_q3_rotation_keeps_mismatch_constant_but_changes_relation() -> None:
    world = simulate_paired_transition_world(
        mode="shared_relational_rotation",
        n_species=5,
        records_per_species=40,
        amplitude=3.0,
        noise_sd=0.0,
        seed=2,
    )
    for sample in world.samples:
        assert np.allclose(pointwise_mismatch(sample), 3.0)
        mismatch_edges = build_mismatch_edges(sample, k=4)
        coupling_edges = build_coupling_edges(sample, k=4)
        assert np.allclose(mismatch_edges.turnover, 0.5)
        assert np.ptp(coupling_edges.turnover) > 0.0


def test_q2_shared_mismatch_front_has_positive_heldout_transfer_smoke() -> None:
    world = simulate_paired_transition_world(
        mode="shared_mismatch_transition",
        n_species=20,
        records_per_species=80,
        amplitude=3.0,
        noise_sd=0.0,
        transition_width=0.15,
        seed=3,
    )
    labels = [sample.species for sample in world.samples]
    train_ids, eval_ids = split_species(labels, eval_fraction=0.5, seed=11)
    by_name = {sample.species: sample for sample in world.samples}

    train_edges = [build_mismatch_edges(by_name[name], k=8) for name in train_ids]
    eval_edges = [build_mismatch_edges(by_name[name], k=8) for name in eval_ids]
    result = transfer_statistic(
        train_edges,
        eval_edges,
        bandwidth=0.35,
        prior_strength=0.0,
        segment_points=5,
    )
    assert result.statistic > 0.0


def test_q1_private_boundaries_are_not_all_at_shared_phase() -> None:
    world = simulate_paired_transition_world(
        mode="private_mismatch_transition",
        n_species=12,
        records_per_species=20,
        seed=4,
    )
    phases = np.array(list(world.boundary_phase.values()), dtype=float)
    assert np.ptp(phases) > 1.0
