from __future__ import annotations

import numpy as np

from ttf.chunked_transfer import prepare_chunked_transfer, score_chunked_batch
from ttf.conditional_transfer import (
    build_conditional_source_sets,
    paired_macro_increment,
    pairwise_geographic_coverage,
    prepare_target_restricted_transfer,
    score_target_restricted_batch,
)
from ttf.core import SpeciesEdges


def _edges(
    name: str,
    shift: float,
    turnover=(0.1, 0.9, 0.3),
) -> SpeciesEdges:
    start = np.array([
        [0.0 + shift, 0.0],
        [1.0 + shift, 0.0],
        [2.0 + shift, 0.0],
    ])
    end = start + np.array([
        [0.5, 0.0],
        [0.5, 0.0],
        [0.5, 0.0],
    ])
    return SpeciesEdges(
        species=name,
        nodes=np.array(
            [[0, 1], [1, 2], [2, 3]],
            dtype=int,
        ),
        start=start,
        end=end,
        midpoint=0.5 * (start + end),
        length=np.linalg.norm(end - start, axis=1),
        turnover=np.asarray(turnover, dtype=float),
    )


def test_geographic_coverage_and_taxonomic_source_sets_are_response_blind() -> None:
    train = [
        _edges("a", 0.0),
        _edges("b", 100.0),
        _edges("c", 0.2),
    ]
    evaluation = [
        _edges("x", 0.1),
        _edges("y", 100.1),
    ]
    coverage = pairwise_geographic_coverage(
        train,
        evaluation,
        radius=2.0,
    )
    assert coverage["x"]["a"] == 1.0
    assert coverage["x"]["b"] == 0.0
    assert coverage["x"]["c"] == 1.0
    assert coverage["y"]["b"] == 1.0

    taxonomy = {
        "a": {"order": "O1"},
        "b": {"order": "O2"},
        "c": {"order": "O1"},
        "x": {"order": "O1"},
        "y": {"order": "O2"},
    }
    frozen = build_conditional_source_sets(
        coverage,
        taxonomy,
        radius=2.0,
        coverage_threshold=0.25,
        min_sources=2,
    )
    assert frozen.radius == 2.0
    assert frozen.geographic["x"] == ("a", "c")
    assert frozen.same_order["x"] == ("a", "c")
    assert frozen.geographic["y"] == ("b",)
    assert frozen.primary_eval_species == ("x",)
    assert frozen.secondary_eval_species == ("x",)


def test_target_restricted_all_sources_matches_baseline_kernel() -> None:
    train = [
        _edges("a", 0.0),
        _edges("b", 0.3, turnover=(0.2, 0.7, 0.4)),
    ]
    evaluation = [
        _edges("x", 0.1, turnover=(0.15, 0.85, 0.25)),
        _edges("y", 0.2, turnover=(0.25, 0.75, 0.35)),
    ]
    width = 3
    train_values = {
        edges.species: np.repeat(
            edges.turnover[:, None],
            width,
            axis=1,
        )
        for edges in train
    }
    eval_values = {
        edges.species: np.repeat(
            edges.turnover[:, None],
            width,
            axis=1,
        )
        for edges in evaluation
    }
    baseline = score_chunked_batch(
        prepare_chunked_transfer(
            train,
            evaluation,
            bandwidth=2.0,
            prior_mean=0.0,
        ),
        train_values,
        eval_values,
    )
    restricted = prepare_target_restricted_transfer(
        train,
        evaluation,
        {
            "x": ("a", "b"),
            "y": ("a", "b"),
        },
        bandwidth=2.0,
        prior_mean=0.0,
    )
    scored = score_target_restricted_batch(
        restricted,
        train_values,
        eval_values,
    )
    assert np.allclose(
        scored.statistics,
        baseline.statistics,
        rtol=0.0,
        atol=1e-14,
    )
    for name in ("x", "y"):
        assert np.allclose(
            scored.species_scores[name],
            baseline.species_scores[name],
            rtol=0.0,
            atol=1e-14,
        )


def test_paired_increment_uses_identical_fixed_targets() -> None:
    train = [
        _edges("a", 0.0),
        _edges("b", 10.0, turnover=(0.9, 0.1, 0.8)),
    ]
    evaluation = [
        _edges("x", 0.1),
        _edges("y", 10.1, turnover=(0.9, 0.1, 0.8)),
    ]
    train_values = {
        edges.species: edges.turnover[:, None]
        for edges in train
    }
    eval_values = {
        edges.species: edges.turnover[:, None]
        for edges in evaluation
    }
    baseline = score_chunked_batch(
        prepare_chunked_transfer(
            train,
            evaluation,
            bandwidth=3.0,
            prior_mean=0.0,
        ),
        train_values,
        eval_values,
    )
    restricted = score_target_restricted_batch(
        prepare_target_restricted_transfer(
            train,
            evaluation,
            {
                "x": ("a",),
                "y": ("b",),
            },
            bandwidth=3.0,
            prior_mean=0.0,
        ),
        train_values,
        eval_values,
    )
    delta = paired_macro_increment(
        restricted,
        baseline,
        ("x", "y"),
    )
    expected = np.mean(
        [
            restricted.species_scores[name]
            - baseline.species_scores[name]
            for name in ("x", "y")
        ],
        axis=0,
    )
    assert np.allclose(delta, expected)
