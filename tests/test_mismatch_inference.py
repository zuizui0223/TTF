import numpy as np

from ttf.core import knn_edges, split_species
from ttf.mismatch import build_coupling_edges, build_mismatch_edges
from ttf.mismatch_inference import paired_heldout_species_bootstrap_test
from ttf.mismatch_simulate import simulate_paired_transition_world
from ttf.transfer import transfer_statistic


def test_joint_paired_inference_matches_two_core_transfer_scores() -> None:
    world = simulate_paired_transition_world(
        mode="shared_mismatch_transition",
        n_species=12,
        records_per_species=30,
        amplitude=2.0,
        noise_sd=0.2,
        transition_width=0.2,
        seed=71,
    )
    labels = [sample.species for sample in world.samples]
    train, evaluation = split_species(labels, eval_fraction=0.5, seed=72)
    by_name = {sample.species: sample for sample in world.samples}

    m_edges = {}
    c_edges = {}
    for name in train + evaluation:
        sample = by_name[name]
        nodes = knn_edges(sample.coordinates, k=4)
        m_edges[name] = build_mismatch_edges(sample, edge_nodes=nodes)
        c_edges[name] = build_coupling_edges(sample, edge_nodes=nodes)

    dense_m = transfer_statistic(
        [m_edges[name] for name in train],
        [m_edges[name] for name in evaluation],
        bandwidth=0.3,
        prior_strength=0.25,
        segment_points=3,
    )
    dense_c = transfer_statistic(
        [c_edges[name] for name in train],
        [c_edges[name] for name in evaluation],
        bandwidth=0.3,
        prior_strength=0.25,
        segment_points=3,
    )

    joint = paired_heldout_species_bootstrap_test(
        world.samples,
        train_species=train,
        eval_species=evaluation,
        k=4,
        bandwidth=0.3,
        prior_strength=0.25,
        segment_points=3,
        n_bootstrap=199,
        seed=73,
        edge_chunk_size=8,
        train_chunk_size=128,
    )

    assert np.isclose(joint.mismatch_statistic, dense_m.statistic, atol=1e-12, rtol=0.0)
    assert np.isclose(joint.coupling_statistic, dense_c.statistic, atol=1e-12, rtol=0.0)
    assert 0.0 < joint.mismatch_bootstrap.p_value <= 1.0
    assert 0.0 < joint.coupling_bootstrap.p_value <= 1.0


def test_joint_inference_keeps_q3_mismatch_null_while_coupling_is_positive() -> None:
    world = simulate_paired_transition_world(
        mode="shared_relational_rotation",
        n_species=16,
        records_per_species=40,
        amplitude=2.0,
        noise_sd=0.0,
        transition_width=0.15,
        seed=81,
    )
    labels = [sample.species for sample in world.samples]
    train, evaluation = split_species(labels, eval_fraction=0.5, seed=82)
    joint = paired_heldout_species_bootstrap_test(
        world.samples,
        train_species=train,
        eval_species=evaluation,
        k=6,
        bandwidth=0.3,
        prior_strength=0.0,
        segment_points=3,
        n_bootstrap=199,
        seed=83,
        edge_chunk_size=8,
        train_chunk_size=256,
    )

    assert np.isclose(joint.mismatch_statistic, 0.0, atol=1e-12, rtol=0.0)
    assert joint.coupling_statistic > 0.0
