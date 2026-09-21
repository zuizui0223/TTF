from __future__ import annotations

from ttf.lepidoptera_host_resource import (
    build_insect_host_footprints,
    jaccard_units,
    select_and_split_species,
    species_list_sha256,
)


def test_host_footprint_union_and_diagnostics():
    pairs=[
        ("A a","h1"),("A a","h2"),("A a","h2"),
        ("B b","h2"),("C c","h3"),
    ]
    units={"h1":["AAA","BBB"],"h2":["BBB","CCC"],"h3":[]}
    footprints,diag=build_insect_host_footprints(pairs,units)
    assert footprints["A a"]==frozenset({"AAA","BBB","CCC"})
    assert footprints["B b"]==frozenset({"BBB","CCC"})
    assert "C c" not in footprints
    assert diag["A a"]["resolved_host_species"]==2
    assert diag["A a"]["hosts_with_primary_native_units"]==2
    assert diag["C c"]["primary_native_wgsrpd3_units"]==0


def test_jaccard_units():
    assert jaccard_units({"A","B","C"},{"B","C","D"})==0.5


def test_panel_split_is_deterministic_disjoint_and_capped():
    species=[f"Species {i:03d}" for i in range(701)]
    selected,train,evaluation=select_and_split_species(species,maximum_species=500)
    selected2,train2,evaluation2=select_and_split_species(reversed(species),maximum_species=500)
    assert selected==selected2
    assert train==train2
    assert evaluation==evaluation2
    assert len(selected)==500 and len(train)==250 and len(evaluation)==250
    assert not (set(train)&set(evaluation))
    assert set(train)|set(evaluation)==set(selected)
    assert len(species_list_sha256(selected))==64
