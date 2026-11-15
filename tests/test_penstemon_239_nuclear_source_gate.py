"""Synthetic public Dryad HAL metadata support, no Penstemon data values."""
import importlib.util
from pathlib import Path
import pytest

P=Path(__file__).resolve().parents[1]/"research"/"prospective"/"penstemon_239_nuclear_source_gate_20261115.py"
spec=importlib.util.spec_from_file_location("penstemon_source",P)
m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)

def test_version_id_from_HAL_self():
    d={"_embedded":{"stash:versions":[
        {"_links":{"self":{"href":"/api/v2/versions/777"}},
         "versionNumber":1,"publicationDate":"2026-11-02"}
    ]}}
    r=m.version_rows(d)
    assert r[0]["id"]==777
    assert r[0]["version"]==1

def test_file_sha256_metadata_and_filename_exact():
    d={"_embedded":{"stash:files":[
        {"path":"2.RAxML.zip","size":650970,
         "digestType":"sha-256","digest":"f"*64,
         "_links":{"self":{"href":"/api/v2/files/12"}}}
    ]}}
    rows=m.file_rows(d)
    assert rows[0]["id"]==12
    assert rows[0]["name"]=="2.RAxML.zip"
    assert rows[0]["sha256"]=="f"*64

def test_no_source_files_should_stop():
    with pytest.raises(ValueError,match="HAL collection"):
        m.file_rows({"_embedded":{"stash:versions":[]}})
