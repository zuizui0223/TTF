import numpy as np

from ttf.external_geographic_opportunity import (
    chord_radius_km,
    directed_occurrence_coverage,
    latlon_to_ecef_km,
)


def test_latlon_to_ecef_preserves_earth_radius():
    xyz = latlon_to_ecef_km(
        np.asarray([0.0, 45.0, -30.0]),
        np.asarray([0.0, 90.0, 120.0]),
    )
    radii = np.sqrt(np.sum(xyz * xyz, axis=1))
    assert np.allclose(radii, 6371.0088)


def test_directed_occurrence_coverage_is_directional():
    source = np.asarray([[0.0, 0.0]])
    target = np.asarray([[0.0, 0.0], [0.0, 10.0]])
    forward = directed_occurrence_coverage(
        target,
        source,
        radius_km=500.0,
    )
    reverse = directed_occurrence_coverage(
        source,
        target,
        radius_km=500.0,
    )
    assert forward == 0.5
    assert reverse == 1.0


def test_chord_radius_matches_arc_threshold_order():
    chord = chord_radius_km(500.0)
    assert 0.0 < chord < 500.0
    near = np.asarray([[0.0, 4.0]])
    far = np.asarray([[0.0, 5.0]])
    source = np.asarray([[0.0, 0.0]])
    assert directed_occurrence_coverage(
        near,
        source,
        radius_km=500.0,
    ) == 1.0
    assert directed_occurrence_coverage(
        far,
        source,
        radius_km=500.0,
    ) == 0.0
