"""Offline tests for source-access versus biological-null distinction."""
import importlib.util
from pathlib import Path
import hashlib
import io
import zipfile
import pytest

PATH=Path(__file__).resolve().parents[1]/"research"/"prospective"/"ttf_crest_costus_original_data_access_gate_20261115.py"
spec=importlib.util.spec_from_file_location("ttf_costus_access",PATH)
m=importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)


def test_checksum_first_csv_header_not_rows():
    raw=b"Species,Pollinator,Count\nA,Bird,12\nB,Bee,3\n"
    item={"sha256":hashlib.sha256(raw).hexdigest(),"size":len(raw),"content_kind":"csv_header","file":"test.csv"}
    r=m.metadata_only(raw,item)
    assert r["header_only"]==["Species","Pollinator","Count"]
    assert r["data_row_values_decoded"] is False
    assert "Bird" not in repr(r) and "Bee" not in repr(r)
    assert "12" not in repr(r)


def test_checksum_drift_stop():
    raw=b"a,b\n1,2\n"
    item={"sha256":"0"*64,"size":len(raw),"content_kind":"csv_header","file":"test.csv"}
    with pytest.raises(ValueError,match="SHA256"):
        m.metadata_only(raw,item)


def test_zip_member_listing_never_opens_archive_data():
    buf=io.BytesIO()
    with zipfile.ZipFile(buf,"w") as z:z.writestr("wild.xlsx",b"unopened bytes")
    raw=buf.getvalue()
    item={"sha256":hashlib.sha256(raw).hexdigest(),"size":len(raw),
          "content_kind":"zip_member_names_only","file":"data.zip"}
    result=m.metadata_only(raw,item)
    assert result["archive_member_names"]==["wild.xlsx"]
    assert result["archive_members_contents_decoded"] is False
    assert "unopened" not in repr(result)
