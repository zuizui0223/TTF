import csv, importlib.util, json
from pathlib import Path

SCRIPT=Path("scripts/aggregate_relational_gbif_precensus.py")

def load():
    spec=importlib.util.spec_from_file_location("relagg",SCRIPT)
    m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m);return m

def test_truthy():
    m=load()
    assert m.truthy("True")
    assert m.truthy("1")
    assert not m.truthy("False")

def test_sha256_path_stable(tmp_path):
    m=load(); p=tmp_path/"x";p.write_text("abc")
    assert m.sha256_path(p)==m.sha256_path(p)
