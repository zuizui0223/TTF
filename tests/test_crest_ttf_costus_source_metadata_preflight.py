"""Synthetic Costus Dryad/OSF source-only metadata shape checks."""
import importlib.util
from pathlib import Path
import pytest

PATH=Path(__file__).resolve().parents[1]/"research"/"prospective"/"ttf_crest_costus_source_metadata_preflight_20261115.py"
spec=importlib.util.spec_from_file_location("ttf_costus_source",PATH)
m=importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)


def test_dryad_hal_version_id_derived_from_link():
    payload={"_embedded":{"stash:versions":[
      {"_links":{"self":{"href":"/api/v2/versions/42"}},
       "versionNumber":2,"publicationDate":"2022-08-05"}
    ]}}
    result=m.parse_versions(payload)
    assert result==[{"id":42,"version_number":2,"publication_date":"2022-08-05",
                     "self":"/api/v2/versions/42"}]


def test_dryad_hal_file_path_and_digest_not_guessed():
    payload={"_embedded":{"stash:files":[
      {"_links":{"self":{"href":"/api/v2/files/99"}},
       "path":"chronogram_v8_2.tre","size":123,
       "digestType":"sha-256","digest":"a"*64}
    ]}}
    result=m.parse_files(payload)
    assert result[0]["filename"]=="chronogram_v8_2.tre"
    assert result[0]["file_id"]==99
    assert result[0]["digest"]=="a"*64
    assert result[0]["bytes"]==123


def test_missing_hal_key_stops_without_fabricating_data():
    with pytest.raises(ValueError,match="HAL collection"):
        m.parse_files({"_embedded":{"stash:versions":[]}})
