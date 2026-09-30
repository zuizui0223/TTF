import runpy
from pathlib import Path


ROOT=Path(__file__).resolve().parents[1]
NS=runpy.run_path(str(ROOT/"scripts/census_palearctic_lgm_subpanel.py"))


RULE={
    "authoritative_subpanel_filter":{
        "taxon":{
            "class":"Insecta",
            "allowed_orders":["Lepidoptera","Hymenoptera","Coleoptera","Hemiptera","Orthoptera","Psocodea"],
            "excluded_coleoptera_families":["Dytiscidae","Elmidae"],
        }
    }
}


def test_strict_terrestrial_core_is_conservative_and_response_independent():
    f=NS["strict_terrestrial"]
    assert f(RULE,class_name="Insecta",order="Lepidoptera",family="Geometridae")
    assert f(RULE,class_name="Insecta",order="Coleoptera",family="Carabidae")
    assert not f(RULE,class_name="Insecta",order="Coleoptera",family="Dytiscidae")
    assert not f(RULE,class_name="Insecta",order="Diptera",family="Drosophilidae")
    assert not f(RULE,class_name="Arachnida",order="Araneae",family="Linyphiidae")


def test_palearctic_fraction_uses_all_selected_localities_as_denominator():
    f=NS["palearctic_fraction"]
    assert f(["Palearctic"]*8+["Nearctic"]*2)==0.8
    assert f(["Palearctic"]*7+["Nearctic"]*2+[""])==0.7
    assert f([])==0.0


def test_role_assignment_is_deterministic_disjoint_and_balanced():
    species=[f"Species {i}" for i in range(41)]
    source,target=NS["assign_roles"](species,"abc123")
    assert len(source)==20
    assert len(target)==21
    assert set(source).isdisjoint(target)
    assert set(source)|set(target)==set(species)
    assert (source,target)==NS["assign_roles"](species,"abc123")
