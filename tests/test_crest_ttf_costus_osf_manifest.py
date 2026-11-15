"""No real Costus OSF data; synthetic OSF JSON metadata typing tests."""
import importlib.util
from pathlib import Path

F=Path(__file__).resolve().parents[1]/"research"/"prospective"/"ttf_crest_costus_osf_manifest_20261115.py"
spec=importlib.util.spec_from_file_location("costus_osf_manifest",F)
m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)


def test_osf_file_metadata_does_not_decode_payload():
    x={"id":"234","attributes":{"name":"scent_data.csv","kind":"file","size":200,
       "extra":{"hashes":{"sha256":"f"*64}}},
       "links":{"download":"https://files.osf.io/download/test"}}
    got=m.meta_file(x)
    assert got["name"]=="scent_data.csv" and got["size"]==200
    assert got["extra_hashes"]["sha256"]=="f"*64
    assert got["download_link"].startswith("https://")
    assert "data" not in got


def test_folder_related_hrefs_accept_string_or_dictionary():
    base={"attributes":{"name":"folder","kind":"folder"}}
    for target in ("https://api.osf.io/v2/folder/",
                   {"href":"https://api.osf.io/v2/folder/"}):
        rec={**base,"relationships":{"files":{"links":{"related":target}}}}
        assert m.meta_file(rec)["related_children"]=="https://api.osf.io/v2/folder/"
