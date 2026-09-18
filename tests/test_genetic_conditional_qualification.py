from __future__ import annotations

import numpy as np

from ttf.genetic_conditional_qualification import (
    frozen_seed,
    make_conditional_order_worlds,
    prepare_conditional_order_qualification,
    score_conditional_order_world_batch,
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
    design = prepare_conditional_order_qualification(
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


def test_frozen_seed_is_stable_and_cell_specific() -> None:
    a = frozen_seed(20260919, "conditional-order", "private_A2", 7)
    b = frozen_seed(20260919, "conditional-order", "private_A2", 7)
    c = frozen_seed(20260919, "conditional-order", "same_order_A2", 7)
    assert a == b
    assert a != c
    assert 0 <= a < 2**64


def test_order_design_uses_same_response_blind_targets_for_paired_increment() -> None:
    design = _fixture()
    assert len(design.eligible_eval_species) == 6
    assert design.eligible_eval_species == design.order_pools.eligible_eval_species
    assert set(design.eligible_eval_species) <= set(design.geographic_pools.eligible_eval_species)
    for target in design.eligible_eval_species:
        assert len(design.order_pools.source_pool[target]) >= 2
        assert set(design.order_pools.source_pool[target]) <= set(
            design.geographic_pools.source_pool[target]
        )
        target_order = design.order_by_species[target]
        assert all(
            design.order_by_species[source] == target_order
            for source in design.order_pools.source_pool[target]
        )


def test_private_and_same_order_worlds_follow_frozen_group_rule() -> None:
    design = _fixture()
    private = make_conditional_order_worlds(
        design,
        cell="private_A2",
        residual_amplitude=2.0,
        absolute_start=0,
        count=1,
        group_mode="private",
    )[0]
    shared = make_conditional_order_worlds(
        design,
        cell="same_order_A2",
        residual_amplitude=2.0,
        absolute_start=0,
        count=1,
        group_mode="same_order",
    )[0]

    assert len(set(private.group_by_species.values())) == len(design.compact.geometries)
    assert set(shared.group_by_species.values()) == {"order:A", "order:B"}
    same_a = [name for name, group in shared.group_by_species.items() if group == "order:A"]
    assert len(same_a) > 1
    assert np.array_equal(
        shared.group_normal["order:A"],
        shared.group_normal[shared.group_by_species[same_a[1]]],
    )


def test_order_world_batch_produces_paired_statistic_and_bootstrap_pvalue() -> None:
    design = _fixture()
    worlds = make_conditional_order_worlds(
        design,
        cell="same_order_A2",
        residual_amplitude=2.0,
        absolute_start=5,
        count=2,
        group_mode="same_order",
    )
    scored = score_conditional_order_world_batch(
        design,
        worlds,
        cell="same_order_A2",
        absolute_start=5,
        bootstrap_resamples=199,
    )
    assert scored.statistics.shape == (2,)
    assert scored.p_values.shape == (2,)
    assert scored.n_eval_species == 6
    assert np.all(np.isfinite(scored.statistics))
    assert np.all((scored.p_values > 0.0) & (scored.p_values <= 1.0))

    matrix = np.vstack([
        scored.species_increments[name]
        for name in design.eligible_eval_species
    ])
    assert np.allclose(scored.statistics, matrix.mean(axis=0), atol=1e-12, rtol=0.0)
    for name in design.eligible_eval_species:
        assert np.allclose(
            scored.species_increments[name],
            scored.order_species_scores[name] - scored.geographic_species_scores[name],
            atol=0.0,
            rtol=0.0,
        )


def test_full_projection_execution_matches_frozen_partial_cache_on_fixture() -> None:
    from ttf.genetic_conditional_qualification import (
        prepare_conditional_order_full_projection_execution,
        score_conditional_order_world_batch_full_projection,
    )

    design = _fixture()
    execution = prepare_conditional_order_full_projection_execution(design)
    worlds = make_conditional_order_worlds(
        design,
        cell="same_order_A2",
        residual_amplitude=2.0,
        absolute_start=17,
        count=2,
        group_mode="same_order",
    )
    expected = score_conditional_order_world_batch(
        design,
        worlds,
        cell="same_order_A2",
        absolute_start=17,
        bootstrap_resamples=199,
    )
    observed = score_conditional_order_world_batch_full_projection(
        design,
        execution,
        worlds,
        cell="same_order_A2",
        absolute_start=17,
        bootstrap_resamples=199,
    )
    assert observed.n_eval_species == expected.n_eval_species
    assert np.allclose(observed.statistics, expected.statistics, atol=1e-12, rtol=0.0)
    assert np.allclose(observed.p_values, expected.p_values, atol=0.0, rtol=0.0)
    for name in design.eligible_eval_species:
        assert np.allclose(
            observed.species_increments[name],
            expected.species_increments[name],
            atol=1e-12,
            rtol=0.0,
        )
        assert np.allclose(
            observed.geographic_species_scores[name],
            expected.geographic_species_scores[name],
            atol=1e-12,
            rtol=0.0,
        )
        assert np.allclose(
            observed.order_species_scores[name],
            expected.order_species_scores[name],
            atol=1e-12,
            rtol=0.0,
        )
