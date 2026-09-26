from __future__ import annotations

from ttf.butterfly_resource_envelope import (
    ResourceEnvelopeDescriptor,
    descriptor_from_sources,
    select_resource_space_pilot,
    species_digest,
    wing_size_proxy,
)


def _descriptor(species: str, host_families: float, host_units: int):
    return ResourceEnvelopeDescriptor(
        species=species,
        host_family_count=host_families,
        host_wgsrpd3_unit_count=host_units,
        wing_size_proxy=20.0,
        voltinism="1",
        canopy_affinity="0",
        edge_affinity="0",
        moisture_affinity="0",
        disturbance_affinity="0",
        resolved_host_species=max(1, int(host_families)),
        hosts_with_primary_native_units=max(1, int(host_families)),
    )


def test_wing_size_proxy_ignores_missing_values():
    row = {
        "WS_L": "10",
        "WS_U": "20",
        "FW_L": "NA",
        "FW_U": "",
    }
    assert wing_size_proxy(row) == 15.0


def test_descriptor_uses_trait_and_host_axes():
    row = {
        "NumberOfHostplantFamilies": "3",
        "Voltinism": "2",
        "CanopyAffinity": "low",
        "EdgeAffinity": "high",
        "MoistureAffinity": "mid",
        "DisturbanceAffinity": "low",
        "WS_L": "10",
        "WS_U": "20",
    }
    d = descriptor_from_sources(
        "Alpha beta",
        row,
        {
            "resolved_host_species": 4,
            "hosts_with_primary_native_units": 3,
            "primary_native_wgsrpd3_units": 17,
        },
    )
    assert d.species == "Alpha beta"
    assert d.host_family_count == 3.0
    assert d.host_wgsrpd3_unit_count == 17
    assert d.wing_size_proxy == 15.0


def test_pilot_selection_is_order_invariant_and_spans_resource_space():
    descriptors = [
        _descriptor("Corner low-low", 1, 1),
        _descriptor("Corner low-high", 1, 100),
        _descriptor("Corner high-low", 100, 1),
        _descriptor("Corner high-high", 100, 100),
        _descriptor("Center a", 10, 10),
        _descriptor("Center b", 12, 8),
    ]
    selected = select_resource_space_pilot(descriptors, pilot_size=4)
    reversed_selected = select_resource_space_pilot(
        list(reversed(descriptors)), pilot_size=4
    )
    assert selected == reversed_selected
    assert set(selected) == {
        "Corner low-low",
        "Corner low-high",
        "Corner high-low",
        "Corner high-high",
    }
    assert len(species_digest(selected)) == 64


def test_pilot_excludes_species_without_native_host_footprint():
    descriptors = [
        _descriptor("A a", 1, 1),
        _descriptor("B b", 2, 2),
        _descriptor("C c", 3, 3),
        _descriptor("No footprint", 100, 0),
    ]
    selected = select_resource_space_pilot(descriptors, pilot_size=3)
    assert "No footprint" not in selected
