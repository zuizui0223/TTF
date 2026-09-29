import importlib.util
from pathlib import Path


SCRIPT=Path("scripts/acquire_palearctic_lgm_occurrences.py")


def load():
    spec=importlib.util.spec_from_file_location("palocc",SCRIPT)
    module=importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


class Realm:
    def intersects(self,point):
        return point.x>=0


def test_realm_filter_applies_before_thinning(monkeypatch):
    m=load()
    calls=[]
    def fake(path,params,retries=8):
        calls.append((path,dict(params)))
        if path=="species/match":
            return {"usageKey":1,"matchType":"EXACT","rank":"SPECIES","canonicalName":"Test species"}
        if params.get("limit")==1:
            return {"count":4}
        return {"results":[
            {"key":1,"decimalLatitude":50,"decimalLongitude":10},
            {"key":2,"decimalLatitude":51,"decimalLongitude":11},
            {"key":3,"decimalLatitude":40,"decimalLongitude":-100},
            {"key":4,"decimalLatitude":41,"decimalLongitude":-101},
        ]}
    monkeypatch.setattr(m,"get_json",fake)
    retained,ledger=m.fetch_species_palearctic(
        "Test species",Realm(),minimum_retained=1,maximum_retained=200
    )
    assert ledger["queried_rows_with_coordinates"]==4
    assert ledger["queried_rows_in_palearctic"]==2
    assert all(float(row["longitude"])>=0 for row in retained)
    assert ledger["status"]=="PASS_OCCURRENCE_GEOMETRY"
