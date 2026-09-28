import numpy as np

from ttf.relational_benchmark import known_truth_design
from ttf.relational_qualification import (
    dyadic_signal_support,
    residualize_primary_relation,
)


def test_known_truth_nested_information_is_recovered_exactly():
    design = known_truth_design(
        n_source=18,
        n_target=16,
        endpoint_retained_fraction=0.35,
        baseline_survival_after_fe=0.55,
        geography_survival_after_baseline=0.70,
        seed=41,
    )

    c1_residual, c1 = residualize_primary_relation(
        design.source,
        design.target,
        design.primary,
        design.baseline_control[:, None],
    )
    c2_residual, c2 = residualize_primary_relation(
        design.source,
        design.target,
        design.primary,
        np.column_stack((
            design.baseline_control,
            design.geographic_control,
        )),
    )

    truth = design.truth
    assert np.isclose(
        c1.source_target_fe_retained_variance_fraction,
        truth["endpoint_retained_fraction"],
        atol=1e-12,
    )
    assert np.isclose(
        c1.control_unique_variance_fraction_after_fe,
        truth["baseline_survival_after_fe"],
        atol=1e-12,
    )
    assert np.isclose(
        c1.total_unique_variance_fraction,
        truth["C1_total_unique_variance_fraction"],
        atol=1e-12,
    )
    assert np.isclose(
        c2.total_unique_variance_fraction,
        truth["C2_total_unique_variance_fraction"],
        atol=1e-12,
    )
    retained = float(np.dot(c2_residual, c2_residual) / np.dot(c1_residual, c1_residual))
    assert np.isclose(
        retained,
        truth["geography_survival_after_baseline"],
        atol=1e-12,
    )


def test_known_truth_concentration_is_detected_without_changing_survival_truth():
    broad = known_truth_design(
        n_source=20,
        n_target=18,
        endpoint_retained_fraction=0.60,
        baseline_survival_after_fe=0.75,
        geography_survival_after_baseline=0.80,
        concentration="broad",
        seed=73,
    )
    concentrated = known_truth_design(
        n_source=20,
        n_target=18,
        endpoint_retained_fraction=0.60,
        baseline_survival_after_fe=0.75,
        geography_survival_after_baseline=0.80,
        concentration="source",
        seed=73,
    )

    def support(design):
        residual, summary = residualize_primary_relation(
            design.source,
            design.target,
            design.primary,
            np.column_stack((
                design.baseline_control,
                design.geographic_control,
            )),
        )
        return summary, dyadic_signal_support(
            design.source,
            design.target,
            residual,
        )

    broad_summary, broad_support = support(broad)
    concentrated_summary, concentrated_support = support(concentrated)

    assert np.isclose(
        broad_summary.total_unique_variance_fraction,
        concentrated_summary.total_unique_variance_fraction,
        atol=1e-12,
    )
    assert np.isclose(
        broad_summary.total_unique_variance_fraction,
        broad.truth["C2_total_unique_variance_fraction"],
        atol=1e-12,
    )
    assert (
        concentrated_support.signal_effective_sources
        < broad_support.signal_effective_sources
    )
    assert (
        concentrated_support.max_source_signal_share
        > broad_support.max_source_signal_share
    )
