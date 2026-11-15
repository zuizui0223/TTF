"""Synthetic original supplementary-archive metadata only; no plant values."""
import importlib.util
from pathlib import Path
from io import BytesIO
import zipfile
import pytest

FILE=Path(__file__).resolve().parents[1]/"research"/"prospective"/"penstemon_europepmc_supplement_archive_gate_20261115.py"
spec=importlib.util.spec_from_file_location("penstemon_europe_supp",FILE)
m=importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)


def test_archive_lists_only_filenames_not_contents():
    buff=BytesIO()
    with zipfile.ZipFile(buff,"w") as z:
        z.writestr("media-1.pdf",b"%PDF-mock accession private_TT_value")
    out=m.zip_support_metadata(buff.getvalue())
    assert out["member_count"]==1
    assert out["member_metadata"][0]["path"]=="media-1.pdf"
    assert "private_TT_value" not in repr(out)
    assert out["no_original_pdf_or_tabular_contents_read"] is True


def test_official_xml_ref_carries_locator_only():
    source=b'''<article xmlns:xlink="http://www.w3.org/1999/xlink"><back>
    <supplementary-material id="supp1" xlink:href="media-1.pdf">
    <table-wrap><tbody><tr><td>CAW-187</td></tr></tbody></table-wrap>
    </supplementary-material></back></article>'''
    out=m.xml_support_metadata(source)
    assert out["root_tag"]=="article"
    assert any(x["refs"].get("href")=="media-1.pdf" for x in out["source_link_records"])
    assert "CAW-187" not in repr(out)
    assert out["no_original_sample_or_phenotypic_values_read"] is True


def test_source_xml_is_rejected_if_not_valid():
    with pytest.raises(Exception):
        m.xml_support_metadata(b"<article><invalid>")
