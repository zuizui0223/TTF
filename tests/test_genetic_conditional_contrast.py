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
