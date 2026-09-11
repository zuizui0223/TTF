import numpy as np

from ttf.heterogeneous_inference import heterogeneous_paired_heldout_species_bootstrap_test
from ttf.heterogeneous_simulate import simulate_q5_world
from ttf.q51_inference import q51_paired_heldout_species_bootstrap_test


def _split():
    train = [f"sp_{i:03d}" for i in range(20)]
    evaluation = [f"sp_{i:03d}" for i in range(20, 40)]
    return train, evaluation


def test_q51_preserves_frozen_ttf_m_and_reports_raw_coupling_audit():
    world = simulate_q5_world(
        response_mode="power_shared_mismatch",
        geometry_profile="matched",
        amplitude=2.0,
        noise_sd=0.4,
        seed=1404,
    )
    train, evaluation = _split()
    kwargs = dict(
        train_species=train,
        eval_species=evaluation,
        graph_fraction=0.15,
        bandwidth=0.2,
        n_bootstrap=99,
        seed=1505,
        edge_chunk_size=32,
        train_chunk_size=2048,
    )
    base = heterogeneous_paired_heldout_species_bootstrap_test(world.samples, **kwargs)
    q51 = q51_paired_heldout_species_bootstrap_test(world.samples, **kwargs)

    assert q51.mismatch_statistic == base.mismatch_statistic
    assert q51.mismatch_species_scores == base.mismatch_species_scores
    assert q51.raw_coupling_statistic == base.coupling_statistic
    assert q51.raw_coupling_species_scores == base.coupling_species_scores
    assert q51.graph_k == base.graph_k
    assert q51.effective_n == base.effective_n
    assert np.isfinite(q51.coupling_statistic)
    assert np.isfinite(q51.coupling_bootstrap.p_value)


def test_q51_coupling_candidate_runs_on_private_relation_shifted_geometry():
    world = simulate_q5_world(
        response_mode="null_private_relation",
        geometry_profile="shifted",
        amplitude=2.0,
        noise_sd=0.0,
        seed=1606,
    )
    train, evaluation = _split()
    result = q51_paired_heldout_species_bootstrap_test(
        world.samples,
        train_species=train,
        eval_species=evaluation,
        graph_fraction=0.15,
        bandwidth=0.2,
        n_bootstrap=99,
        seed=1707,
        edge_chunk_size=32,
        train_chunk_size=2048,
    )
    assert np.isfinite(result.coupling_statistic)
    assert len(result.coupling_species_scores) == 20
    assert 0.0 <= result.coupling_bootstrap.p_value <= 1.0
