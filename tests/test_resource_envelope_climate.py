from __future__ import annotations

import numpy as np
import pytest

from ttf.resource_envelope_climate import (
    climate_mismatch,
    deterministic_interior_points,
    observed_unit_split,
    probability_greater,
)


def test_observed_unit_split_is_deterministic_disjoint_and_complete():
    units = ["AAA", "BBB", "CCC", "DDD", "EEE", "FFF"]
    train, evaluation = observed_unit_split("Species alpha", units)
    train2, evaluation2 = observed_unit_split("Species alpha", reversed(units))
    assert train == train2
    assert evaluation == evaluation2
    assert set(train).isdisjoint(evaluation)
    assert set(train) | set(evaluation) == set(units)
    assert train and evaluation


def test_climate_mismatch_is_zero_at_training_center():
    train = np.asarray(
        [
            [0.0, 0.0],
            [2.0, 4.0],
            [1.0, 2.0],
        ]
    )
    units = np.asarray([[1.0, 2.0], [2.0, 4.0]])
    mismatch = climate_mismatch(units, train)
    assert np.isclose(mismatch[0], 0.0)
    assert mismatch[1] > mismatch[0]


def test_probability_greater_handles_ties():
    assert probability_greater([2.0, 3.0], [1.0, 2.0]) == 0.875
    assert probability_greater([], [1.0]) is None


def test_deterministic_interior_points_cover_multipolygon_components():
    shapely = pytest.importorskip("shapely")
    from shapely.geometry import MultiPolygon, Point, Polygon

    a = Polygon([(0, 0), (1, 0), (1, 1), (0, 1)])
    b = Polygon([(10, 10), (11, 10), (11, 11), (10, 11)])
    geom = MultiPolygon([a, b])
    points = deterministic_interior_points(geom, maximum_points=8)
    assert points == deterministic_interior_points(geom, maximum_points=8)
    assert 2 <= len(points) <= 8
    assert all(geom.covers(Point(x, y)) for x, y in points)
    assert any(x < 2 for x, _ in points)
    assert any(x > 9 for x, _ in points)
