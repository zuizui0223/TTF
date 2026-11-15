"""Synthetic OSF nested folder graph contract; no biological files."""
import importlib.util
from pathlib import Path

FILE=Path(__file__).resolve().parents[1]/"research"/"prospective"/"ttf_crest_costus_osf_nested_manifest_20261115.py"
spec=importlib.util.spec_from_file_location("costus_osf_nested",FILE)
m=importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)


def test_osf_link_can_be_string_or_object():
    for val in ("https://api.osf.io/v2/folder/xyz/",
                {"href":"https://api.osf.io/v2/folder/xyz/"}):
        r={"relationships":{"files":{"links":{"related":val}}}}
        assert m.folder_href(r)=="https://api.osf.io/v2/folder/xyz/"


def test_documented_root_ids_are_distinct():
    r=m.FOLDERS
    assert r["2025_flower_scent"]["node"]=="2ap4k"
    assert r["2025_phylogenies"]["node"]=="2ap4k"
    assert r["2022_phylogenomics_pollination"]["node"]=="fkj2x"
    assert len({(v["node"],v["folder"]) for v in r.values()})==3
