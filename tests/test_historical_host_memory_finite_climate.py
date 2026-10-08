import numpy as np
import pytest

from scripts.census_historical_host_memory_finite_climate import census_finite_host_support


def test_support_count_is_unique_union_of_host_native_units():
    candidates = [
        {"species": "Alpha beta", "accepted_host_ids": "400000.0;900"},
        {"species": "Gamma delta", "accepted_host_ids": "123"},
    ]
    native = [
        {"accepted_plant_name_id": "4e+05", "area_code_l3": "AA"},
        {"accepted_plant_name_id": "900", "area_code_l3": "AA"},
        {"accepted_plant_name_id": "900", "area_code_l3": "BB"},
        {"accepted_plant_name_id": "123", "area_code_l3": "CC"},
    ]
    support = [{"area_code_l3": s} for s in ("AA", "AA", "BB", "BB", "CC")]
    result = census_finite_host_support(
        candidates, native, support, np.array([True, False, True, True, False])
    )
    assert [x["host_climate_points"] for x in result] == [3, 0]


def test_invalid_support_mask_refuses_silent_alignment():
    with pytest.raises(RuntimeError, match="length mismatch"):
        census_finite_host_support(
            [{"species": "A b", "accepted_host_ids": "1"}],
            [{"accepted_plant_name_id": "1", "area_code_l3": "A"}],
            [{"area_code_l3": "A"}],
            np.array([True, False]),
        )


def test_any_empty_native_unit_or_support_area_abstains():
    cand = [{"species": "A b", "accepted_host_ids": "1"}]
    with pytest.raises(RuntimeError, match="empty native"):
        census_finite_host_support(
            cand, [{"accepted_plant_name_id": "1", "area_code_l3": ""}],
            [{"area_code_l3": "A"}], np.array([True]))
    with pytest.raises(RuntimeError, match="empty WGSRPD3"):
        census_finite_host_support(
            cand, [{"accepted_plant_name_id": "1", "area_code_l3": "A"}],
            [{"area_code_l3": ""}], np.array([True]))
