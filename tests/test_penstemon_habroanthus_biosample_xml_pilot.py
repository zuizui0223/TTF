"""Synthetic Penstemon BioSample XML identity pilot tests (offline only)."""
import importlib.util
from pathlib import Path
import xml.etree.ElementTree as ET
import pytest

PATH=Path(__file__).resolve().parents[1]/"research"/"prospective"/"penstemon_habroanthus_biosample_xml_pilot_20261115.py"
spec=importlib.util.spec_from_file_location("penstemon_2026_biosample_xml",PATH)
mod=importlib.util.module_from_spec(spec)
spec.loader.exec_module(mod)


def test_minimal_xml_extracts_taxon_and_only_attribute_names():
    xml=b"""<SAMPLE_SET><SAMPLE accession="SAMN999" alias="CAW123">
        <TITLE>Some original plant sample</TITLE>
        <SAMPLE_NAME><TAXON_ID>12345</TAXON_ID>
        <SCIENTIFIC_NAME>Penstemon barbatus</SCIENTIFIC_NAME></SAMPLE_NAME>
        <SAMPLE_ATTRIBUTES>
           <SAMPLE_ATTRIBUTE><TAG>collector</TAG><VALUE>Not for reuse</VALUE></SAMPLE_ATTRIBUTE>
           <SAMPLE_ATTRIBUTE><TAG>voucher</TAG><VALUE>Voucher-123</VALUE></SAMPLE_ATTRIBUTE>
        </SAMPLE_ATTRIBUTES>
    </SAMPLE></SAMPLE_SET>"""
    x=mod.inspect_xml(xml)
    assert x["source_xml_record_count"]==1
    assert x["scientific_name"]=="Penstemon barbatus"
    assert x["taxon_id"]=="12345"
    assert x["submitter_alias_present"] is True
    assert x["title_present"] is True
    assert x["available_attribute_field_names"]==["collector","voucher"]
    assert "Voucher-123" not in repr(x)
    assert "Not for reuse" not in repr(x)
    assert x["no_sequence_or_trait_values_read"] is True


def test_bare_sample_and_missing_taxon_is_unknown_not_species_absence():
    x=mod.inspect_xml(b"<SAMPLE accession='SAMN1'><SAMPLE_NAME/></SAMPLE>")
    assert x["source_xml_record_count"]==1
    assert x["scientific_name"] is None
    assert x["taxon_id"] is None


def test_multi_record_response_refused():
    with pytest.raises(ValueError,match="record count"):
        mod.inspect_xml(b"<SAMPLE_SET><SAMPLE/><SAMPLE/></SAMPLE_SET>")


def test_malformed_xml_refused():
    with pytest.raises(ET.ParseError):
        mod.inspect_xml(b"<SAMPLE><")


def test_deterministic_five_sample_cap():
    assert mod.LIMIT==5
    assert mod.STUDY=="PRJNA1479764"
