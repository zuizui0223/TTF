import json
from pathlib import Path

import numpy as np

from ttf.relational_qualification import (
    detectability_surface,
    dyadic_signal_support,
    minimum_detectable_effects,
    relation_repeatability,
    residualize_primary_relation,
)


ROOT = Path(__file__).resolve().parents[1]


def complete_bipartite(ns=12, nt=10):
    source = np.repeat(np.arange(ns), nt)
    target = np.tile(np.arange(nt), ns)
    return source, target


def test_relation_repeatability_distinguishes_stable_from_noisy_relations():
    rng = np.random.default_rng(4)
    truth = rng.normal(size=80)
    stable = truth[None, :] + rng.normal(
        scale=0.05,
        size=(12, len(truth)),
    )
    noisy = truth[None, :] + rng.normal(
        scale=2.0,
        size=(12, len(truth)),
    )
    stable_result = relation_repeatability(stable)
    noisy_result = relation_repeatability(noisy)
    assert stable_result.repeatability_icc > 0.95
    assert noisy_result.repeatability_icc < stable_result.repeatability_icc
    assert (
        stable_result.median_within_dyad_sd
        < noisy_result.median_within_dyad_sd
    )


def test_nonredundancy_reports_relation_variance_not_explained_by_controls():
    source, target = complete_bipartite()
    rng = np.random.default_rng(8)
    control = rng.normal(size=len(source))
    almost_redundant = control + rng.normal(
        scale=0.03,
        size=len(source),
    )
    residual, summary = residualize_primary_relation(
        source,
        target,
        (almost_redundant - almost_redundant.mean())
        / almost_redundant.std(),
        control[:, None],
    )
    assert summary.control_unique_variance_fraction_after_fe < 0.01
    assert summary.vif_like_after_fe > 100
    assert np.std(residual) < 0.2


def test_information_loss_decomposition_is_exact():
    source, target = complete_bipartite(11, 9)
    rng = np.random.default_rng(81)
    primary = rng.normal(size=len(source))
    controls = np.column_stack((
        0.4 * primary + rng.normal(size=len(source)),
        rng.normal(size=len(source)),
    ))
    _, summary = residualize_primary_relation(
        source,
        target,
        (primary - primary.mean()) / primary.std(),
        controls,
    )
    assert np.isclose(
        summary.total_unique_variance_fraction,
        summary.source_target_fe_retained_variance_fraction
        * summary.control_unique_variance_fraction_after_fe,
        rtol=1e-12,
        atol=1e-12,
    )


def test_dyadic_signal_support_detects_endpoint_concentration():
    source, target = complete_bipartite(10, 8)
    broad = np.ones(len(source))
    broad[::2] = -1
    concentrated = np.where(source == 0, 10.0, 0.1)
    broad_support = dyadic_signal_support(
        source,
        target,
        broad,
    )
    concentrated_support = dyadic_signal_support(
        source,
        target,
        concentrated,
    )
    assert broad_support.signal_effective_sources > 9.5
    assert concentrated_support.signal_effective_sources < 2.0
    assert concentrated_support.max_source_signal_share > 0.9


def test_detectability_surface_is_continuous_characterization_not_gate():
    source, target = complete_bipartite()
    rng = np.random.default_rng(12)
    primary = rng.normal(size=len(source))
    control = rng.normal(size=len(source))
    predictors = np.column_stack((primary, control))
    surface = detectability_surface(
        source,
        target,
        predictors,
        effects=[0.0, 0.8],
        private_amplitudes=[0.0],
        alpha=0.05,
        worlds_per_cell=200,
        source_intercept_sd=0.15,
        target_intercept_sd=0.15,
        dyad_noise_sd=0.25,
        fixed_coefficients=[0.0, 0.1],
        response_transform="identity",
        seed_namespace="test-ttf-q",
        block_size=50,
    )
    by_effect = {
        float(row["effect"]): row
        for row in surface
    }
    assert by_effect[0.0]["rejection_rate"] < 0.15
    assert by_effect[0.8]["rejection_rate"] > 0.95
    mde = minimum_detectable_effects(
        surface,
        target_power=0.8,
    )
    assert mde["0"] == 0.8


def test_ttf_q_contract_does_not_reopen_closed_v03():
    contract = json.loads(
        (ROOT / "docs/supporting/ttf_q_v0.1.json").read_text()
    )
    final_state = json.loads(
        (
            ROOT
            / "benchmarks/frozen/relational_program_final_state_v0.3.json"
        ).read_text()
    )
    assert final_state["final_state"] == "CLOSED_NO_EVALUABLE_TEST"
    assert (
        contract["boundary"]["closed_v03_final_state"]
        == final_state["final_state"]
    )
    assert contract["boundary"]["v03_reopening_authorized"] is False
    assert (
        contract["boundary"]["v03_slot_replacement_authorized"]
        is False
    )
    assert (
        contract["boundary"]["genetic_response_opening_authorized"]
        is False
    )
    assert (
        contract["interpretation_policy"]["binary_pass_fail_forbidden"]
        is True
    )
