from __future__ import annotations
import importlib.util
from pathlib import Path

SCRIPT=Path("scripts/census_lepidoptera_trait_transfer_panel.py")
spec=importlib.util.spec_from_file_location("lep_census",SCRIPT)
mod=importlib.util.module_from_spec(spec); assert spec.loader is not None; spec.loader.exec_module(mod)

def rows(order,n=12,support=5):
    return [{"order":order,"locality_index":str(i),"min_endpoint_disjoint_training_edges":str(support)} for i in range(n)]

def test_lepidoptera_metadata_filter():
    x={"good":rows("Lepidoptera"),"few":rows("Lepidoptera",11),"weak":rows("Lepidoptera",12,4),"fly":rows("Diptera")}
    assert mod.eligible_leps(x)==("good",)

def test_trait_funnel_is_staged_not_all_or_nothing():
    x={"WS_U":"4.2","Voltinism":"M","NumberOfHostplantFamilies":"2","CanopyAffinity":"Closed canopy","EdgeAffinity":"NA","MoistureAffinity":"Mesic","DisturbanceAffinity":"Associated"}
    f=mod.flags(x)
    assert f["wing_size"] and f["voltinism"] and f["host_breadth"] and f["habitat_any"]
    assert not f["habitat_complete4"] and not f["major_complete"]

def test_deterministic_split_is_disjoint_exhaustive():
    x=tuple(f"sp{i}" for i in range(11))
    a,b=mod.split(x,"tag"); c,d=mod.split(tuple(reversed(x)),"tag")
    assert (a,b)==(c,d)
    assert set(a).isdisjoint(b) and set(a)|set(b)==set(x)
    assert (len(a),len(b))==(5,6)

def test_schema_is_locked():
    assert mod.SCHEMA=="ttf_lepidoptera_trait_transfer_census_v0.1"
