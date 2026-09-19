from __future__ import annotations

import numpy as np

from ttf.chunked_transfer import prepare_chunked_transfer, score_chunked_batch
from ttf.conditional_transfer import (
    infer_conditioning_increment,
    prepare_target_source_pools,
    score_conditioning_increment_batch,
    score_pool_difference_batch,
    score_target_conditioned_batch,
)
from ttf.core import SpeciesEdges


def _edge(name: str, offset: float, turnover: np.ndarray) -> SpeciesEdges:
    x = np.arange(8, dtype=float) + offset
    coordinates = np.column_stack([x, np.zeros_like(x)])
    nodes = np.column_stack([np.arange(7), np.arange(1, 8)])
    start, end = coordinates[nodes[:, 0]], coordinates[nodes[:, 1]]
    return SpeciesEdges(
        species=name,
        nodes=nodes,
        start=start,
        end=end,
        midpoint=0.5 * (start + end),
        length=np.linalg.norm(end - start, axis=1),
        turnover=np.asarray(turnover, dtype=float),
    )


def _fixture():
    pattern = np.linspace(0.0, 1.0, 7)
    train = [
        _edge("near_a", 0.0, pattern),
        _edge("near_b", 0.2, pattern[::-1]),
        _edge("far_a", 100.0, pattern),
        _edge("far_b", 200.0, pattern[::-1]),
        _edge("far_c", 300.0, pattern),
    ]
    evaluation = [
        _edge(f"eval_{i}", 0.1 + 0.01 * i, pattern) for i in range(6)
    ]
    prepared = prepare_chunked_transfer(
        train,
        evaluation,
        bandwidth=2.0,
        prior_strength=0.25,
        prior_mean=0.0,
        segment_points=3,
    )
    return train, evaluation, prepared


def test_source_pools_are_geometry_only_and_freeze_support() -> None:
    train, evaluation, prepared = _fixture()
    pools = prepare_target_source_pools(
        prepared,
        {x.species: x.midpoint for x in train},
        {x.species: x.midpoint for x in evaluation},
        support_radius=5.0,
        minimum_target_coverage=0.25,
        minimum_source_species=2,
    )
    assert pools.eligible_eval_species == tuple(x.species for x in evaluation)
    assert all(
        pools.source_pool[x.species] == ("near_a", "near_b") for x in evaluation
    )
    assert all(
        pools.pair_coverage[x.species]["far_a"] == 0.0 for x in evaluation
    )


def test_same_group_pool_requires_response_blind_label_support() -> None:
    train, evaluation, prepared = _fixture()
    train_group = {
        x.species: ("A" if x.species.startswith("near") else "B") for x in train
    }
    eval_group = {x.species: "A" for x in evaluation}
    pools = prepare_target_source_pools(
        prepared,
        {x.species: x.midpoint for x in train},
        {x.species: x.midpoint for x in evaluation},
        support_radius=5.0,
        minimum_target_coverage=0.25,
        minimum_source_species=2,
        train_group=train_group,
        eval_group=eval_group,
        require_same_group=True,
    )
    assert all(
        pools.source_pool[x.species] == ("near_a", "near_b") for x in evaluation
    )


def test_conditioned_batch_matches_baseline_when_every_source_is_selected() -> None:
    train, evaluation, prepared = _fixture()
    pools = prepare_target_source_pools(
        prepared,
        {x.species: x.midpoint for x in train},
        {x.species: x.midpoint for x in evaluation},
        support_radius=1000.0,
        minimum_target_coverage=0.0,
        minimum_source_species=5,
    )
    train_y = {x.species: x.turnover[:, None] for x in train}
    eval_y = {x.species: x.turnover[:, None] for x in evaluation}
    baseline = score_chunked_batch(
        prepared, train_y, eval_y, edge_chunk_size=3, train_chunk_size=5
    )
    conditioned = score_target_conditioned_batch(
        prepared, pools, train_y, eval_y, edge_chunk_size=3, train_chunk_size=5
    )
    assert np.allclose(conditioned.statistics, baseline.statistics, atol=1e-12, rtol=0)
    for name in pools.eligible_eval_species:
        assert np.allclose(
            conditioned.species_scores[name],
            baseline.species_scores[name],
            atol=1e-12,
            rtol=0,
        )


def test_increment_is_paired_by_target_and_bootstraps_species() -> None:
    train, evaluation, prepared = _fixture()
    pools = prepare_target_source_pools(
        prepared,
        {x.species: x.midpoint for x in train},
        {x.species: x.midpoint for x in evaluation},
        support_radius=5.0,
        minimum_target_coverage=0.25,
        minimum_source_species=2,
    )
    train_y = {x.species: x.turnover[:, None] for x in train}
    eval_y = {x.species: x.turnover[:, None] for x in evaluation}
    scored = score_conditioning_increment_batch(prepared, pools, train_y, eval_y)
    assert scored.statistics.shape == (1,)
    assert len(scored.species_increments) == 6
    for name in pools.eligible_eval_species:
        expected = (
            scored.conditional_species_scores[name]
            - scored.baseline_species_scores[name]
        )
        assert np.allclose(scored.species_increments[name], expected, atol=0, rtol=0)
    inferred = infer_conditioning_increment(scored, n_bootstrap=199, seed=11)
    assert np.isclose(
        inferred.statistic,
        np.mean(list(inferred.species_increments.values())),
    )
    assert 0.0 < inferred.bootstrap.p_value <= 1.0


def test_order_like_pool_difference_is_paired_on_supported_targets() -> None:
    train, evaluation, prepared = _fixture()
    train_mid = {x.species: x.midpoint for x in train}
    eval_mid = {x.species: x.midpoint for x in evaluation}
    geographic = prepare_target_source_pools(
        prepared,
        train_mid,
        eval_mid,
        support_radius=5.0,
        minimum_target_coverage=0.25,
        minimum_source_species=2,
    )
    train_group = {
        x.species: ("A" if x.species.startswith("near") else "B") for x in train
    }
    eval_group = {x.species: "A" for x in evaluation}
    order_like = prepare_target_source_pools(
        prepared,
        train_mid,
        eval_mid,
        support_radius=5.0,
        minimum_target_coverage=0.25,
        minimum_source_species=2,
        train_group=train_group,
        eval_group=eval_group,
        require_same_group=True,
    )
    train_y = {x.species: x.turnover[:, None] for x in train}
    eval_y = {x.species: x.turnover[:, None] for x in evaluation}
    scored = score_pool_difference_batch(
        prepared, geographic, order_like, train_y, eval_y
    )
    assert scored.eligible_eval_species == order_like.eligible_eval_species
    assert scored.statistics.shape == (1,)
    for name in scored.eligible_eval_species:
        assert np.allclose(
            scored.species_increments[name],
            scored.conditional_species_scores[name]
            - scored.baseline_species_scores[name],
        )


def test_cached_conditioned_scoring_matches_uncached_exactly() -> None:
    from ttf.cached_chunked_transfer import prepare_cached_chunked_transfer
    from ttf.conditional_transfer import (
        prepare_cached_target_conditioned_transfer,
        score_cached_conditioning_increment_batch,
        score_cached_target_conditioned_batch,
    )

    train, evaluation, prepared = _fixture()
    pools = prepare_target_source_pools(
        prepared,
        {x.species: x.midpoint for x in train},
        {x.species: x.midpoint for x in evaluation},
        support_radius=5.0,
        minimum_target_coverage=0.25,
        minimum_source_species=2,
    )
    train_y = {
        x.species: np.column_stack([x.turnover, x.turnover[::-1]])
        for x in train
    }
    eval_y = {
        x.species: np.column_stack([x.turnover, x.turnover[::-1]])
        for x in evaluation
    }
    uncached = score_target_conditioned_batch(prepared, pools, train_y, eval_y)
    cached_cond = prepare_cached_target_conditioned_transfer(prepared, pools)
    cached = score_cached_target_conditioned_batch(cached_cond, train_y, eval_y)
    assert np.allclose(cached.statistics, uncached.statistics, atol=1e-12, rtol=0)
    for name in pools.eligible_eval_species:
        assert np.allclose(
            cached.species_scores[name],
            uncached.species_scores[name],
            atol=1e-12,
            rtol=0,
        )

    baseline_cache = prepare_cached_chunked_transfer(prepared)
    increment = score_cached_conditioning_increment_batch(
        baseline_cache, cached_cond, train_y, eval_y
    )
    assert increment.statistics.shape == (2,)


def test_conditioned_pool_preserves_total_training_species_weight_mass() -> None:
    from ttf.conditional_transfer import prepare_cached_target_conditioned_transfer

    train, evaluation, prepared = _fixture()
    pools = prepare_target_source_pools(
        prepared,
        {x.species: x.midpoint for x in train},
        {x.species: x.midpoint for x in evaluation},
        support_radius=5.0,
        minimum_target_coverage=0.50,
        minimum_source_species=2,
    )
    cached = prepare_cached_target_conditioned_transfer(prepared, pools)

    full_mass = float(np.sum(prepared.train_weights))
    assert np.isclose(full_mass, float(len(prepared.train_species)))
    for target in pools.eligible_eval_species:
        active = cached.active_indices[target]
        scaled_mass = float(
            np.sum(prepared.train_weights[active]) * cached.pool_weight_scale[target]
        )
        assert np.isclose(scaled_mass, full_mass, atol=1e-12, rtol=0.0)
        assert np.isclose(
            cached.pool_weight_scale[target],
            len(prepared.train_species) / len(pools.source_pool[target]),
            atol=0.0,
            rtol=0.0,
        )


def test_full_projection_cache_matches_existing_cached_scorer_to_frozen_tolerance() -> None:
    from ttf.conditional_transfer import (
        prepare_cached_target_conditioned_transfer,
        prepare_fully_cached_target_conditioned_transfer,
        score_cached_target_conditioned_batch,
        score_fully_cached_target_conditioned_batch,
    )

    train, evaluation, prepared = _fixture()
    pools = prepare_target_source_pools(
        prepared,
        {x.species: x.midpoint for x in train},
        {x.species: x.midpoint for x in evaluation},
        support_radius=5.0,
        minimum_target_coverage=0.50,
        minimum_source_species=2,
    )
    partial = prepare_cached_target_conditioned_transfer(
        prepared, pools, edge_chunk_size=3, train_chunk_size=5
    )
    full = prepare_fully_cached_target_conditioned_transfer(partial)
    train_y = {
        x.species: np.column_stack([x.turnover, x.turnover[::-1]])
        for x in train
    }
    eval_y = {
        x.species: np.column_stack([x.turnover, x.turnover[::-1]])
        for x in evaluation
    }
    expected = score_cached_target_conditioned_batch(partial, train_y, eval_y)
    observed = score_fully_cached_target_conditioned_batch(full, train_y, eval_y)

    assert np.allclose(observed.statistics, expected.statistics, atol=1e-12, rtol=0.0)
    assert np.array_equal(observed.n_eval_species, expected.n_eval_species)
    for name in pools.eligible_eval_species:
        assert np.allclose(
            observed.species_scores[name],
            expected.species_scores[name],
            atol=1e-12,
            rtol=0.0,
        )


def test_full_projection_cache_can_materialize_predeclared_target_subset() -> None:
    from ttf.conditional_transfer import (
        prepare_cached_target_conditioned_transfer,
        prepare_fully_cached_target_conditioned_transfer,
        score_fully_cached_target_conditioned_batch,
    )

    train, evaluation, prepared = _fixture()
    pools = prepare_target_source_pools(
        prepared,
        {x.species: x.midpoint for x in train},
        {x.species: x.midpoint for x in evaluation},
        support_radius=5.0,
        minimum_target_coverage=0.50,
        minimum_source_species=2,
    )
    partial = prepare_cached_target_conditioned_transfer(prepared, pools)
    subset = tuple(pools.eligible_eval_species[:3])
    full = prepare_fully_cached_target_conditioned_transfer(
        partial, eval_species=subset
    )
    train_y = {x.species: x.turnover[:, None] for x in train}
    eval_y = {x.species: x.turnover[:, None] for x in evaluation}
    scored = score_fully_cached_target_conditioned_batch(full, train_y, eval_y)
    assert full.eval_species == subset
    assert set(scored.species_scores) == set(subset)
    assert np.array_equal(scored.n_eval_species, np.full(1, len(subset)))


def test_different_group_pool_is_disjoint_from_same_group_pool() -> None:
    train, evaluation, prepared = _fixture()
    train_mid = {x.species: x.midpoint for x in train}
    eval_mid = {x.species: x.midpoint for x in evaluation}
    train_group = {
        "near_a": "A",
        "near_b": "A",
        "far_a": "B",
        "far_b": "B",
        "far_c": "B",
    }
    eval_group = {x.species: "A" for x in evaluation}
    same = prepare_target_source_pools(
        prepared,
        train_mid,
        eval_mid,
        support_radius=1000.0,
        minimum_target_coverage=0.0,
        minimum_source_species=2,
        train_group=train_group,
        eval_group=eval_group,
        require_same_group=True,
    )
    different = prepare_target_source_pools(
        prepared,
        train_mid,
        eval_mid,
        support_radius=1000.0,
        minimum_target_coverage=0.0,
        minimum_source_species=2,
        train_group=train_group,
        eval_group=eval_group,
        require_different_group=True,
    )
    assert same.require_same_group is True
    assert same.require_different_group is False
    assert different.require_same_group is False
    assert different.require_different_group is True
    for target in same.eligible_eval_species:
        assert same.source_pool[target] == ("near_a", "near_b")
        assert different.source_pool[target] == ("far_a", "far_b", "far_c")
        assert set(same.source_pool[target]).isdisjoint(different.source_pool[target])


def test_source_pool_rejects_conflicting_group_rules() -> None:
    train, evaluation, prepared = _fixture()
    train_group = {x.species: "A" for x in train}
    eval_group = {x.species: "A" for x in evaluation}
    import pytest
    with pytest.raises(ValueError, match="both same and different"):
        prepare_target_source_pools(
            prepared,
            {x.species: x.midpoint for x in train},
            {x.species: x.midpoint for x in evaluation},
            train_group=train_group,
            eval_group=eval_group,
            require_same_group=True,
            require_different_group=True,
        )
