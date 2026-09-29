from scripts.fetch_wwf_palearctic_realm import canonical_query_url, validate_geojson


def test_wwf_realm_query_is_frozen_to_palearctic_wgs84_geojson():
    url=canonical_query_url()
    assert "realmcode%3D%27PA%27" in url
    assert "outSR=4326" in url
    assert "f=geojson" in url


def test_wwf_realm_fetch_validation_rejects_non_palearctic():
    import json
    good={
        "type":"FeatureCollection",
        "features":[{
            "type":"Feature",
            "properties":{"realm":"Palearctic","realmcode":"PA"},
            "geometry":{"type":"Polygon","coordinates":[[[0,0],[1,0],[1,1],[0,1],[0,0]]]},
        }],
    }
    validate_geojson(json.dumps(good).encode())
    bad={
        "type":"FeatureCollection",
        "features":[{
            "type":"Feature",
            "properties":{"realm":"Nearctic","realmcode":"NA"},
            "geometry":{"type":"Polygon","coordinates":[[[0,0],[1,0],[1,1],[0,1],[0,0]]]},
        }],
    }
    import pytest
    with pytest.raises(RuntimeError,match="non-Palearctic"):
        validate_geojson(json.dumps(bad).encode())
