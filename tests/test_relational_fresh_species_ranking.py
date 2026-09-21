import importlib.util
from pathlib import Path

SCRIPT=Path("scripts/freeze_relational_fresh_species_ranking.py")

def load():
    spec=importlib.util.spec_from_file_location("relrank",SCRIPT)
    m=importlib.util.module_from_spec(spec); spec.loader.exec_module(m); return m

def test_rank_is_deterministic_and_source_bound():
    m=load()
    assert m.rank_key("abc","Species a")==m.rank_key("abc","Species a")
    assert m.rank_key("abc","Species a")!=m.rank_key("def","Species a")

def test_tag_is_fresh_relational_namespace():
    m=load()
    assert m.TAG=="relational-ttf-v0.1"
