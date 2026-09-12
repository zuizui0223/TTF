import numpy as np

from ttf.mismatch import PairedSpeciesSample
from ttf.q51_inference import q51_paired_heldout_species_bootstrap_test
from ttf.training_leverage_diagnostic import paired_training_leverage_diagnostic


def _small_samples(response_variant: int = 0):
    samples = []
    base_theta = np.linspace(0.0, 2.0 * np.pi, 12, endpoint=False)
    for i in range(12):
        theta = np.mod(base_theta + 0.025 * i, 2.0 * np.pi)
        coordinates = np.column_stack((np.cos(theta), np.sin(theta)))
        phase = 0.35 * i
        if response_variant == 0:
            state_a = 0.15 * np.cos(theta + 0.1 * i)
            state_b = np.tanh(np.sin(theta - phase) / 0.25) + 0.05 * np.cos(2.0 * theta + i)
        else:
            state_a = 0.4 * np.sin(3.0 * theta + i)
            state_b = 1.7 * np.cos(theta + 0.2 * i) + 0.2 * np.sin(4.0 * theta - phase)
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


def _kwargs():
    train, evaluation = _split()
    return dict(
        train_species=train,
        eval_species=evaluation,
        graph_fraction=0.20,
        bandwidth=0.45,
        prior_strength=0.25,
        segment_points=5,
    )


def test_training_leverage_full_scores_reproduce_frozen_q5_q51_estimands():
    samples = _small_samples()
    kwargs = _kwargs()
    q51 = q51_paired_heldout_species_bootstrap_test(
        samples,
        **kwargs,
        n_bootstrap=99,
        seed=4401,
        edge_chunk_size=16,
        train_chunk_size=256,
    )
    diagnostic = paired_training_leverage_diagnostic(samples, **kwargs)

    assert np.isclose(
        diagnostic.mismatch.raw_statistics[0], q51.mismatch_statistic, atol=1e-12, rtol=0.0
    )
    assert np.isclose(
        diagnostic.coupling.raw_statistics[0], q51.coupling_statistic, atol=1e-12, rtol=0.0
    )
    assert diagnostic.graph_k == q51.graph_k
    assert diagnostic.effective_n == q51.effective_n


def test_training_system_concentration_is_response_blind_on_fixed_geometry():
    first = paired_training_leverage_diagnostic(_small_samples(0), **_kwargs())
    second = paired_training_leverage_diagnostic(_small_samples(1), **_kwargs())

    assert np.isclose(
        first.mean_effective_training_system_count,
        second.mean_effective_training_system_count,
        atol=1e-12,
        rtol=0.0,
    )
    assert np.isclose(
        first.mean_dominant_training_system_share,
        second.mean_dominant_training_system_share,
        atol=1e-12,
        rtol=0.0,
    )


def test_training_leverage_reports_all_exact_leaveouts_and_finite_endpoints():
    result = paired_training_leverage_diagnostic(_small_samples(), **_kwargs())
    train, _ = _split()

    assert 1.0 <= result.mean_effective_training_system_count <= len(train)
    assert 1.0 / len(train) <= result.mean_dominant_training_system_share <= 1.0
    for estimand in (result.mismatch, result.coupling):
        assert set(estimand.loo_statistics_by_training_system) == set(train)
        for name in train:
            assert np.isfinite(estimand.loo_statistics_by_training_system[name]).all()
        assert np.isfinite(estimand.raw_statistics).all()
        assert np.isfinite(estimand.loo_score_sd).all()
        assert np.isfinite(estimand.loo_score_mean).all()
        assert np.isfinite(estimand.loo_score_min).all()
        assert np.isfinite(estimand.mean_abs_loo_score_shift).all()
        assert np.isfinite(estimand.mean_prediction_loo_sd).all()
        assert np.all(estimand.loo_score_sd >= 0.0)
        assert np.all(estimand.mean_abs_loo_score_shift >= 0.0)
        assert np.all(estimand.mean_prediction_loo_sd >= 0.0)
