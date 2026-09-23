import importlib.util
from pathlib import Path

SCRIPT=Path("scripts/freeze_relational_gbif_precensus.py")

def load():
    spec=importlib.util.spec_from_file_location("relgbif",SCRIPT)
    m=importlib.util.module_from_spec(spec); spec.loader.exec_module(m); return m

def test_inspect_species_count_gate(monkeypatch):
    m=load()
    def fake(path, params, retries=4):
        if path=="species/match":
            return {"usageKey": 123, "matchType":"EXACT", "canonicalName":"Alpha beta", "rank":"SPECIES"}
        if path=="occurrence/search":
            return {"count": 42, "results":[]}
        raise AssertionError(path)
    monkeypatch.setattr(m,"get_json",fake)
    row=m.inspect_species("Alpha beta",start_year=2010,end_year=2026,minimum_count=30,allow_canonical_fuzzy=False)
    assert row["match_accepted"] is True
    assert row["occurrence_count_2010_endyear"]==42
    assert row["count_ge_minimum"] is True
    assert row["status"]=="COUNT_GE_MINIMUM"

def test_rejected_match_never_queries_occurrence(monkeypatch):
    m=load()
    calls=[]
    def fake(path, params, retries=4):
        calls.append(path)
        return {"usageKey": 123, "matchType":"FUZZY", "canonicalName":"Other species", "rank":"SPECIES"}
    monkeypatch.setattr(m,"get_json",fake)
    row=m.inspect_species("Alpha beta",start_year=2010,end_year=2026,minimum_count=30,allow_canonical_fuzzy=False)
    assert row["status"]=="REJECTED_GBIF_TAXON_MATCH"
    assert calls==["species/match"]
