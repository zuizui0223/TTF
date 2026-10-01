from __future__ import annotations

import numpy as np
import pytest

from ttf.codistribution import (
    coverage_fraction,
    geographic_coopportunity,
    geographic_coopportunity_table,
)


def test_coverage_fraction_is_directional_for_nested_support():
    small = np.array([[0.0, 0.0], [1.0, 0.0]])
    large = np.array([[0.0, 0.0], [1.0, 0.0], [20.0, 0.0], [30.0, 0.0]])

    assert coverage_fraction(small, large, radius=2.0) == 1.0
    assert coverage_fraction(large, small, radius=2.0) == 0.5


def test_primary_coopportunity_is_symmetric_minimum():
    source = np.array([[0.0, 0.0], [1.0, 0.0]])
    target = np.array([[0.0, 0.0], [1.0, 0.0], [20.0, 0.0], [30.0, 0.0]])

    forward = geographic_coopportunity("s", "t", source, target, radius=2.0)
    reverse = geographic_coopportunity("t", "s", target, source, radius=2.0)

    assert forward.source_to_target_coverage == 1.0
    assert forward.target_to_source_coverage == 0.5
    assert forward.symmetric_min_coverage == 0.5
    assert reverse.symmetric_min_coverage == forward.symmetric_min_coverage


def test_distant_lineages_have_zero_coopportunity():
    source = np.array([[0.0, 0.0], [1.0, 0.0]])
    target = np.array([[100.0, 0.0], [101.0, 0.0]])

    result = geographic_coopportunity("s", "t", source, target, radius=10.0)

    assert result.source_to_target_coverage == 0.0
    assert result.target_to_source_coverage == 0.0
    assert result.symmetric_min_coverage == 0.0


def test_identical_support_has_complete_coopportunity():
    points = np.array([[0.0, 0.0], [3.0, 4.0], [10.0, 1.0]])

    result = geographic_coopportunity("s", "t", points, points.copy(), radius=0.01)

    assert result.symmetric_min_coverage == 1.0
    assert result.centroid_distance == 0.0


def test_table_retains_all_directed_dyads_without_pool_selection():
    source = {
        "a": np.array([[0.0, 0.0], [1.0, 0.0]]),
        "b": np.array([[10.0, 0.0], [11.0, 0.0]]),
    }
    target = {
        "x": np.array([[0.0, 0.0], [1.0, 0.0]]),
        "y": np.array([[10.0, 0.0], [11.0, 0.0]]),
        "z": np.array([[50.0, 0.0], [51.0, 0.0]]),
    }

    rows = geographic_coopportunity_table(source, target, radius=2.0)

    assert len(rows) == 6
    assert {(row.source, row.target) for row in rows} == {
        ("a", "x"),
        ("b", "x"),
        ("a", "y"),
        ("b", "y"),
        ("a", "z"),
        ("b", "z"),
    }


def test_table_excludes_self_when_source_and_target_sets_overlap():
    points = {
        "a": np.array([[0.0, 0.0], [1.0, 0.0]]),
        "b": np.array([[2.0, 0.0], [3.0, 0.0]]),
    }

    rows = geographic_coopportunity_table(points, points, radius=5.0)

    assert {(row.source, row.target) for row in rows} == {("a", "b"), ("b", "a")}


def test_invalid_geometry_is_rejected():
    with pytest.raises(ValueError):
        coverage_fraction(np.array([0.0, 1.0]), np.array([[0.0, 1.0]]))

    with pytest.raises(ValueError):
        coverage_fraction(
            np.array([[0.0, 1.0]]),
            np.array([[0.0, 1.0, 2.0]]),
        )

    with pytest.raises(ValueError):
        geographic_coopportunity(
            "s",
            "t",
            np.array([[0.0, np.nan]]),
            np.array([[0.0, 0.0]]),
        )
