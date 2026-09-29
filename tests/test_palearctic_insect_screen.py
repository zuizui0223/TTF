import json

import pytest

from ttf.palearctic_insect_screen import (
    apply_gbif_admissibility,
    deterministic_role_split,
    geometry_contains_or_boundary,
    select_realm_geometries,
    summarize_species_realm,
)


def feature_collection():
    return {
        "type": "FeatureCollection",
        "features": [
            {
                "type": "Feature",
                "properties": {"realm": "Palearctic", "realmcode": "PA"},
                "geometry": {
                    "type": "Polygon",
                    "coordinates": [[
                        [0, 0], [10, 0], [10, 10], [0, 10], [0, 0]
                    ]],
                },
            },
            {
                "type": "Feature",
                "properties": {"realm": "Nearctic", "realmcode": "NA"},
                "geometry": {
                    "type": "Polygon",
                    "coordinates": [[
                        [20, 0], [30, 0], [30, 10], [20, 10], [20, 0]
                    ]],
                },
            },
        ],
    }


def test_boundary_and_hole_semantics():
    polygon = {
        "type": "Polygon",
        "coordinates": [
            [[0,0],[10,0],[10,10],[0,10],[0,0]],
            [[4,4],[6,4],[6,6],[4,6],[4,4]],
        ],
    }
    assert geometry_contains_or_boundary(polygon, lon=0, lat=5)
    assert geometry_contains_or_boundary(polygon, lon=3, lat=3)
    assert not geometry_contains_or_boundary(polygon, lon=5, lat=5)
    assert geometry_contains_or_boundary(polygon, lon=4, lat=5)


def test_realm_screen_separates_realm_and_terrestrial_eligibility():
    geoms = select_realm_geometries(feature_collection())
    rows = []
    for species, order, family, points in [
        ("Terra alpha", "Lepidoptera", "Noctuidae", [(1,1),(2,2),(3,3),(4,4),(25,5)]),
        ("Water beta", "Coleoptera", "Dytiscidae", [(1,2),(2,3),(3,4),(4,5),(5,6)]),
        ("Mixed gamma", "Hymenoptera", "Formicidae", [(1,1),(2,2),(25,2),(26,3),(27,4)]),
    ]:
        for lon, lat in points:
            rows.append({
                "species": species,
                "class": "Insecta",
                "order": order,
                "family": family,
                "longitude": lon,
                "latitude": lat,
            })
    out = {
        row.species: row
        for row in summarize_species_realm(
            rows,
            geometries=geoms,
            allowed_orders=["Lepidoptera","Hymenoptera","Coleoptera"],
            excluded_coleoptera_families=["Dytiscidae"],
            minimum_fraction=0.80,
            minimum_localities=4,
        )
    }
    assert out["Terra alpha"].palearctic_fraction == pytest.approx(0.80)
    assert out["Terra alpha"].biologically_eligible
    assert out["Water beta"].realm_eligible
    assert not out["Water beta"].terrestrial_core
    assert not out["Mixed gamma"].realm_eligible


def test_gbif_filter_precedes_deterministic_role_split():
    geoms = select_realm_geometries(feature_collection())
    rows = []
    for index, species in enumerate(["A a","B b","C c","D d","E e","F f"]):
        for j in range(4):
            rows.append({
                "species": species,
                "class": "Insecta",
                "order": "Lepidoptera",
                "family": "Noctuidae",
                "longitude": 1 + j,
                "latitude": 1 + index * 0.1,
            })
    summaries = summarize_species_realm(
        rows,
        geometries=geoms,
        allowed_orders=["Lepidoptera"],
        excluded_coleoptera_families=[],
        minimum_fraction=0.80,
        minimum_localities=4,
    )
    ledger={"A a":30,"B b":29,"C c":40,"D d":50,"E e":31,"F f":0}
    eligible=apply_gbif_admissibility(summaries,ledger,minimum_occurrences=30)
    assert eligible == ("A a","C c","D d","E e")
    source,target=deterministic_role_split(eligible,namespace="test")
    assert len(source)==2 and len(target)==2
    assert set(source).isdisjoint(target)
    assert set(source)|set(target)==set(eligible)


def test_missing_gbif_species_fails_closed():
    geoms=select_realm_geometries(feature_collection())
    summaries=summarize_species_realm(
        [{
            "species":"A a","class":"Insecta","order":"Lepidoptera","family":"Noctuidae",
            "longitude":x,"latitude":x,
        } for x in (1,2,3,4)],
        geometries=geoms,
        allowed_orders=["Lepidoptera"],
        excluded_coleoptera_families=[],
        minimum_fraction=0.8,
        minimum_localities=4,
    )
    with pytest.raises(ValueError, match="missing biologically eligible species"):
        apply_gbif_admissibility(summaries,{},minimum_occurrences=30)
