from __future__ import annotations

import json
from pathlib import Path

import numpy as np

from ttf.genetic_conditional_contrast import (
    BOOTSTRAP_SEED_TAG,
    WORLD_SEED_TAG,
    frozen_contrast_seed,
    make_conditional_order_contrast_worlds,
    prepare_conditional_order_contrast,
    prepare_conditional_order_contrast_execution,
    score_conditional_order_contrast_world_batch,
    summarize_conditional_order_contrast_breadth,
)
from ttf.genetic_geometry import prepare_density_scaled_genetic_geometry


def _fixture():
    theta = np.linspace(0.0, 2.0 * np.pi, 16, endpoint=False)
    base = np.column_stack([np.cos(theta), np.sin(theta)])
    geometries = {}
    for i in range(12):
        shift = np.array([0.02 * (i // 4), 0.02 * (i // 4)])
        geometries[f"sp{i:02d}"] = prepare_density_scaled_genetic_geometry(base + shift)
    names = tuple(sorted(geometries))
    train = names[:6]
    evaluation = names[6:]
    order = {
        name: ("A" if int(name[-2:]) % 2 == 0 else "B")
        for name in names
    }
    design = prepare_conditional_order_contrast(
        geometries,
        order,
        train_species=train,
        eval_species=evaluation,
        bandwidth=2.0,
        support_radius=5.0,
        minimum_target_coverage=0.25,
        minimum_source_species=2,
        edge_chunk_size=8,
        train_chunk_size=128,
    )
    return design


def test_order_contrast_uses_disjoint_geography_matched_source_pools() -> None:
    design = _fixture()
    assert len(design.eligible_eval_species) == 6
    for target in design.eligible_eval_species:
        same = design.same_order_pools.source_pool[target]
        different = design.different_order_pools.source_pool[target]
        assert len(same) >= 2
        assert len(different) >= 2
        assert set(same).isdisjoint(different)
        target_order = design.order_by_species[target]
        assert all(design.order_by_species[source] == target_order for source in same)
        assert all(design.order_by_species[source] != target_order for source in different)


def test_order_contrast_seed_namespace_is_separate_from_v01() -> None:
    assert WORLD_SEED_TAG == "conditional-order-contrast-v02"
    assert BOOTSTRAP_SEED_TAG == "conditional-order-contrast-v02-bootstrap"
    a = frozen_contrast_seed(20260919, WORLD_SEED_TAG, "same_order_A2", 0)
    b = frozen_contrast_seed(20260919, "conditional-order", "same_order_A2", 0)
    assert a != b


def test_order_contrast_world_batch_is_paired_same_minus_different() -> None:
    design = _fixture()
    execution = prepare_conditional_order_contrast_execution(design)
    worlds = make_conditional_order_contrast_worlds(
        design,
        cell="same_order_A2",
        residual_amplitude=2.0,
        absolute_start=3,
        count=2,
        group_mode="same_order",
    )
    scored = score_conditional_order_contrast_world_batch(
        design,
        execution,
        worlds,
        cell="same_order_A2",
        absolute_start=3,
        bootstrap_resamples=199,
    )
    assert scored.statistics.shape == (2,)
    assert scored.p_values.shape == (2,)
    assert scored.n_eval_species == 6
    matrix = np.vstack([
        scored.species_increments[name]
        for name in design.eligible_eval_species
    ])
    assert np.allclose(scored.statistics, matrix.mean(axis=0), atol=1e-12, rtol=0.0)
    for name in design.eligible_eval_species:
        assert np.allclose(
            scored.species_increments[name],
            scored.same_order_species_scores[name]
            - scored.different_order_species_scores[name],
            atol=0.0,
            rtol=0.0,
        )


def test_v02_protocol_freezes_disjoint_formal_seed_namespace_and_closed_firewall() -> None:
    rule = json.loads(
        Path("docs/supporting/genetic_conditional_order_contrast_qualification_rule_v0.2.json").read_text()
    )
    assert rule["status"] == "FROZEN_BEFORE_V02_FORMAL_QUALIFICATION_AND_ANY_CONDITIONAL_EMPIRICAL_IDENTITY"
    assert rule["frozen_estimator"]["minimum_same_order_source_species"] == 5
    assert rule["frozen_estimator"]["minimum_different_order_source_species"] == 5
    assert rule["development_panel"]["jointly_supported_eval_species"] == 52
    assert rule["confirmatory_reserve"]["jointly_supported_eval_species_before_character_masks"] == 77
    assert rule["synthetic_worlds"]["world_seed_tag"] == WORLD_SEED_TAG
    assert rule["synthetic_worlds"]["bootstrap_seed_tag"] == BOOTSTRAP_SEED_TAG
    assert all(value is False for value in rule["outcome_firewall"].values())
    assert rule["development_pilot_disclosure"]["status"] == "METHOD_DEVELOPMENT_ONLY_NOT_FORMAL_QUALIFICATION"
    assert rule["execution"]["world_batch_size"] == 20
    assert "absolute replicate index" in rule["execution"]["sharding_semantics"]


def test_v02_breadth_diagnostic_keeps_species_equal_primary_unchanged() -> None:
    design = _fixture()
    values = {
        name: np.asarray([float(i), float(i + 1)], dtype=float)
        for i, name in enumerate(design.eligible_eval_species)
    }
    diagnostic = summarize_conditional_order_contrast_breadth(
        design,
        values,
        minimum_targets_per_order=2,
    )
    matrix = np.vstack([values[name] for name in design.eligible_eval_species])
    assert np.allclose(
        diagnostic.species_equal_statistics,
        matrix.mean(axis=0),
        atol=0.0,
        rtol=0.0,
    )
    assert diagnostic.order_species_counts == {"A": 3, "B": 3}
    assert diagnostic.included_orders == ("A", "B")
    assert diagnostic.excluded_orders == ()
    assert diagnostic.top_order == "A"
    assert diagnostic.top_order_fraction == 0.5
    for label in ("A", "B"):
        retained = [
            values[name]
            for name in design.eligible_eval_species
            if design.order_by_species[name] != label
        ]
        assert np.allclose(
            diagnostic.leave_one_order_out_statistics[label],
            np.mean(np.vstack(retained), axis=0),
            atol=0.0,
            rtol=0.0,
        )


def test_v02_breadth_rule_is_response_blind_and_descriptive_only() -> None:
    rule = json.loads(
        Path(
            "docs/supporting/genetic_conditional_order_contrast_breadth_rule_v0.2.json"
        ).read_text()
    )
    assert rule["status"].startswith("FROZEN_BEFORE_FORMAL_V02")
    assert rule["breadth_diagnostic"]["inference"] == "descriptive_only_no_p_value"
    dev = rule["response_blind_support_census"]["development"]
    conf = rule["response_blind_support_census"]["confirmatory"]
    assert dev["jointly_supported_targets"] == 52
    assert dev["order_counts"]["Lepidoptera"] == 27
    assert dev["largest_order_fraction"] == 27 / 52
    assert conf["jointly_supported_targets_before_character_masks"] == 77
    assert conf["order_counts"]["Lepidoptera"] == 46
    assert conf["largest_order_fraction"] == 46 / 77
    assert all(value is False for value in rule["outcome_firewall"].values())


def test_survivor_requalification_failure_keeps_identity_closed() -> None:
    receipt = json.loads(
        Path(
            "benchmarks/frozen/genetic_conditional_order_contrast_survivor_requalification_v0.2.json"
        ).read_text()
    )
    assert receipt["status"] == "NOT_EVALUABLE_CONDITIONAL_ORDER_CONTRAST_V02_SURVIVOR_GEOMETRY"
    assert receipt["character_mask_stage"]["survivors"] == 214
    assert receipt["survivor_support"]["jointly_supported_eval_species"] == 63
    assert receipt["frozen_gate"]["maximum_private_rejections_compatible_with_type1_pass"] == 36
    assert receipt["deterministic_partial_execution"]["cells"]["private_A1"]["rejections"] == 45
    assert receipt["deterministic_partial_execution"]["cells"]["private_A2"]["rejections"] == 55
    assert receipt["deterministic_partial_execution"]["cells"]["private_A3"]["rejections"] == 44
    assert receipt["decision"]["type1_gate"] == "FAIL_IRREVERSIBLY_BEFORE_500"
    assert receipt["decision"]["confirmatory_identity_opening_authorized"] is False
    assert all(value is False for value in receipt["outcome_firewall"].values())
