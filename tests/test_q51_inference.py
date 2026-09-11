import numpy as np

from ttf.heterogeneous_inference import heterogeneous_paired_heldout_species_bootstrap_test
from ttf.mismatch import PairedSpeciesSample
from ttf.q51_inference import q51_paired_heldout_species_bootstrap_test


def _small_samples():
    samples = []
    base_theta = np.linspace(0.0, 2.0 * np.pi, 12, endpoint=False)
    for i in range(12):
        theta = np.mod(base_theta + 0.025 * i, 2.0 * np.pi)
        coordinates = np.column_stack((np.cos(theta), np.sin(theta)))
        phase = 0.35 * i
        state_a = np.zeros(len(theta), dtype=float)
        state_b = np.tanh(np.sin(theta - phase) / 0.25) + 0.05 * np.cos(2.0 * theta + i)
        samples.append(
            PairedSpeciesSample(
                species=f"sp_{i:03d}",
                coordinates=coordinates,
                state_a=state_a,
                state_b=state_b,
            )
        )
    return tuple(samples)


def _split():
    train = [f"sp_{i:03d}" for i in range(6)]
    evaluation = [f"sp_{i:03d}" for i in range(6, 12)]
    return train, evaluation


def test_q51_preserves_frozen_ttf_m_and_reports_raw_coupling_audit():
    samples = _small_samples()
    train, evaluation = _split()
    kwargs = dict(
        train_species=train,
        eval_species=evaluation,
        graph_fraction=0.20,
        bandwidth=0.45,
        n_bootstrap=99,
        seed=1505,
        edge_chunk_size=16,
        train_chunk_size=256,
    )
    base = heterogeneous_paired_heldout_species_bootstrap_test(samples, **kwargs)
    q51 = q51_paired_heldout_species_bootstrap_test(samples, **kwargs)

    assert q51.mismatch_statistic == base.mismatch_statistic
    assert q51.mismatch_species_scores == base.mismatch_species_scores
    assert q51.raw_coupling_statistic == base.coupling_statistic
    assert q51.raw_coupling_species_scores == base.coupling_species_scores
    assert q51.graph_k == base.graph_k
    assert q51.effective_n == base.effective_n
    assert np.isfinite(q51.coupling_statistic)
    assert np.isfinite(q51.coupling_bootstrap.p_value)


def test_q51_coupling_candidate_runs_with_system_disjoint_split():
    samples = _small_samples()
    train, evaluation = _split()
    result = q51_paired_heldout_species_bootstrap_test(
        samples,
        train_species=train,
        eval_species=evaluation,
        graph_fraction=0.20,
        bandwidth=0.45,
        n_bootstrap=99,
        seed=1707,
        edge_chunk_size=16,
        train_chunk_size=256,
    )
    assert np.isfinite(result.coupling_statistic)
    assert len(result.coupling_species_scores) == 6
    assert 0.0 <= result.coupling_bootstrap.p_value <= 1.0
