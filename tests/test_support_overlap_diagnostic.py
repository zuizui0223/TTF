import numpy as np

from ttf.mismatch import PairedSpeciesSample
from ttf.q51_inference import q51_paired_heldout_species_bootstrap_test
from ttf.support_overlap_diagnostic import paired_support_overlap_diagnostic


def _small_samples(response_variant: int = 0):
    samples = []
    base_theta = np.linspace(0.0, 2.0 * np.pi, 12, endpoint=False)
    for i in range(12):
        theta = np.mod(base_theta + 0.025 * i, 2.0 * np.pi)
        coordinates = np.column_stack((np.cos(theta), np.sin(theta)))
        phase = 0.35 * i
        if response_variant == 0:
            state_a = 0.15 * np.cos(theta + 0.1 * i)
            state_b = (
                np.tanh(np.sin(theta - phase) / 0.25)
                + 0.05 * np.cos(2.0 * theta + i)
            )
        else:
            # Deliberately different responses on exactly the same geometry.
            state_a = 0.4 * np.sin(3.0 * theta + i)
            state_b = (
                1.7 * np.cos(theta + 0.2 * i)
                + 0.2 * np.sin(4.0 * theta - phase)
            )
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
        edge_chunk_size=16,
        train_chunk_size=256,
    )


def test_support_diagnostic_raw_scores_reproduce_frozen_q5_q51_estimands():
    samples = _small_samples()
    kwargs = _kwargs()
    q51 = q51_paired_heldout_species_bootstrap_test(
        samples,
        **kwargs,
        n_bootstrap=99,
        seed=3301,
    )
    diagnostic = paired_support_overlap_diagnostic(samples, **kwargs)

    assert np.isclose(
        diagnostic.mismatch.raw_statistics[0],
        q51.mismatch_statistic,
        atol=1e-12,
        rtol=0.0,
    )
    assert np.isclose(
        diagnostic.coupling.raw_statistics[0],
        q51.coupling_statistic,
        atol=1e-12,
        rtol=0.0,
    )
    assert diagnostic.graph_k == q51.graph_k
    assert diagnostic.effective_n == q51.effective_n


def test_support_summaries_are_response_blind_on_fixed_geometry():
    kwargs = _kwargs()
    first = paired_support_overlap_diagnostic(_small_samples(0), **kwargs)
    second = paired_support_overlap_diagnostic(_small_samples(1), **kwargs)

    for left, right in ((first.mismatch, second.mismatch), (first.coupling, second.coupling)):
        assert np.isclose(left.mean_opportunity, right.mean_opportunity, atol=1e-12, rtol=0.0)
        assert np.isclose(
            left.mean_prior_fraction,
            right.mean_prior_fraction,
            atol=1e-12,
            rtol=0.0,
        )
        assert np.isclose(
            left.low_support_fraction,
            right.low_support_fraction,
            atol=1e-12,
            rtol=0.0,
        )
        assert left.species_mean_opportunity == right.species_mean_opportunity
        assert left.species_mean_prior_fraction == right.species_mean_prior_fraction
        assert left.species_low_support_fraction == right.species_low_support_fraction


def test_support_diagnostic_keeps_m_and_c_on_identical_geometry_and_finite_endpoints():
    result = paired_support_overlap_diagnostic(_small_samples(), **_kwargs())

    assert np.isclose(
        result.mismatch.mean_opportunity,
        result.coupling.mean_opportunity,
        atol=1e-12,
        rtol=0.0,
    )
    assert np.isclose(
        result.mismatch.mean_prior_fraction,
        result.coupling.mean_prior_fraction,
        atol=1e-12,
        rtol=0.0,
    )
    assert np.isclose(
        result.mismatch.low_support_fraction,
        result.coupling.low_support_fraction,
        atol=1e-12,
        rtol=0.0,
    )

    for estimand in (result.mismatch, result.coupling):
        assert np.isfinite(estimand.raw_statistics).all()
        assert np.isfinite(estimand.opportunity_partial_statistics).all()
        assert np.isfinite(estimand.alignment_statistics).all()
        assert 0.0 <= estimand.mean_prior_fraction <= 1.0
        assert 0.0 <= estimand.low_support_fraction <= 1.0
        assert estimand.mean_opportunity >= 0.0
