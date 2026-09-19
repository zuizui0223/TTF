from __future__ import annotations

import numpy as np

from ttf.genetic_conditional_geometry_match import (
    DEV_EVALUATION_SEED_TAG,
    DEV_REFERENCE_SEED_TAG,
    calibrate_conditional_private_envelope,
    frozen_v03_development_seed,
    make_conditional_geometry_matched_order_worlds,
    prepare_conditional_geometry_matched_order,
    prepare_conditional_geometry_matched_order_execution,
    score_conditional_geometry_matched_order_world_batch,
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
    return prepare_conditional_geometry_matched_order(
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


def test_v03_geometry_matching_uses_equal_disjoint_source_counts() -> None:
    design = _fixture()
    assert len(design.eligible_eval_species) == 6
    for target in design.eligible_eval_species:
        same = design.matched_pools.same_group_pools.source_pool[target]
        different = design.matched_pools.different_group_pools.source_pool[target]
        assert len(same) == len(different) >= 2
        assert set(same).isdisjoint(different)
        target_order = design.order_by_species[target]
        assert all(design.order_by_species[name] == target_order for name in same)
        assert all(design.order_by_species[name] != target_order for name in different)


def test_v03_development_seed_namespaces_are_disjoint_and_locked() -> None:
    a = frozen_v03_development_seed(
        20260919, DEV_REFERENCE_SEED_TAG, "private_A2", 0
    )
    b = frozen_v03_development_seed(
        20260919, DEV_EVALUATION_SEED_TAG, "private_A2", 0
    )
    assert a != b
    assert a == frozen_v03_development_seed(
        20260919, DEV_REFERENCE_SEED_TAG, "private_A2", 0
    )
    import pytest
    with pytest.raises(ValueError, match="not authorized"):
        frozen_v03_development_seed(20260919, "future-formal-tag", "private_A2", 0)


def test_v03_world_score_is_same_order_minus_geometry_matched_different_order() -> None:
    design = _fixture()
    execution = prepare_conditional_geometry_matched_order_execution(design)
    worlds = make_conditional_geometry_matched_order_worlds(
        design,
        cell="same_order_A2",
        residual_amplitude=2.0,
        absolute_start=4,
        count=2,
        group_mode="same_order",
    )
    scored = score_conditional_geometry_matched_order_world_batch(
        design, execution, worlds
    )
    assert scored.statistics.shape == (2,)
    assert scored.n_eval_species == 6
    matrix = np.vstack(
        [scored.species_increments[name] for name in design.eligible_eval_species]
    )
    assert np.allclose(scored.statistics, matrix.mean(axis=0), atol=1e-12, rtol=0.0)
    for name in design.eligible_eval_species:
        assert np.allclose(
            scored.species_increments[name],
            scored.same_order_species_scores[name]
            - scored.different_order_species_scores[name],
            atol=0.0,
            rtol=0.0,
        )


def test_v03_private_envelope_uses_worst_case_independent_reference() -> None:
    references = {
        "private_A0p5": np.linspace(-0.20, 0.20, 99),
        "private_A1": np.linspace(-0.10, 0.30, 99),
        "private_A2": np.linspace(-0.05, 0.35, 99),
        "private_A3": np.linspace(0.00, 0.40, 99),
    }
    observed = np.asarray([0.10, 0.25], dtype=float)
    calibrated = calibrate_conditional_private_envelope(observed, references)
    assert calibrated.statistics.shape == (2,)
    assert calibrated.p_values.shape == (2,)
    assert np.all((calibrated.p_values > 0.0) & (calibrated.p_values <= 1.0))
    # The least-favourable component must be reported rather than silently
    # selecting the most favourable private-null amplitude.
    assert calibrated.least_favourable_private_cell == (
        "private_A3",
        "private_A3",
    )
