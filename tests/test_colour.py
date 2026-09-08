from __future__ import annotations

import numpy as np
import pytest

from ttf.colour import four_component_colour_jsd


def test_colour_jsd_is_zero_on_identity_and_symmetric() -> None:
    p = np.array([0.1, 0.2, 0.3, 0.4])
    q = np.array([0.7, 0.1, 0.1, 0.1])
    assert four_component_colour_jsd(p, p) == pytest.approx(0.0, abs=1e-15)
    assert four_component_colour_jsd(p, q) == pytest.approx(four_component_colour_jsd(q, p), abs=1e-15)


def test_colour_jsd_handles_zero_components_and_known_disjoint_limit() -> None:
    p = np.array([1.0, 0.0, 0.0, 0.0])
    q = np.array([0.0, 1.0, 0.0, 0.0])
    assert four_component_colour_jsd(p, q) == pytest.approx(np.log(2.0), rel=0.0, abs=1e-15)


def test_colour_jsd_rejects_invalid_vectors_without_renormalizing() -> None:
    good = np.array([0.25, 0.25, 0.25, 0.25])
    with pytest.raises(ValueError):
        four_component_colour_jsd(np.array([0.2, 0.2, 0.2, 0.2]), good)
    with pytest.raises(ValueError):
        four_component_colour_jsd(np.array([-0.1, 0.4, 0.3, 0.4]), good)
    with pytest.raises(ValueError):
        four_component_colour_jsd(np.array([0.5, 0.5, 0.0]), good)


def test_jsd_and_sqrt_jsd_have_identical_strict_order() -> None:
    anchor = np.array([1.0, 0.0, 0.0, 0.0])
    others = [
        np.array([0.9, 0.1, 0.0, 0.0]),
        np.array([0.6, 0.4, 0.0, 0.0]),
        np.array([0.0, 1.0, 0.0, 0.0]),
    ]
    jsd = np.array([four_component_colour_jsd(anchor, q) for q in others])
    root = np.sqrt(jsd)
    assert np.array_equal(np.argsort(jsd), np.argsort(root))
