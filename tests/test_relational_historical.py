import numpy as np
import pytest

from ttf.relational_historical import (
    historical_delta,
    logical_asset_basename,
    panel_order,
    panel_roles,
)


def test_frozen_historical_asset_names():
    assert logical_asset_basename("bio01", -190) == "CHELSA_TraCE21K_bio01_-190_V1.0.tif"
    assert logical_asset_basename("bio15", 20) == "CHELSA_TraCE21K_bio15_20_V1.0.tif"
    with pytest.raises(ValueError):
        logical_asset_basename("bio99", 20)


def test_historical_delta_is_lgm_minus_present():
    past = np.asarray([[2.0, 3.0, 4.0, 5.0]])
    present = np.asarray([[1.0, 1.5, 2.0, 2.5]])
    assert np.allclose(historical_delta(past, present), [[1.0, 1.5, 2.0, 2.5]])


def test_history_panel_hashes_are_deterministic_and_disjoint():
    names = [f"Species {i}" for i in range(11)]
    order = panel_order(names, "a" * 64)
    assert order == panel_order(names, "a" * 64)
    source, target = panel_roles(order, "a" * 64)
    assert not (set(source) & set(target))
    assert set(source) | set(target) == set(names)
