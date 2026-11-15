"""Synthetic Penstemon ENA BioSample source-field and metadata grain checks."""
import importlib.util
from pathlib import Path
import pytest

P=Path(__file__).resolve().parents[1]/"research"/"prospective"/"penstemon_habroanthus_ena_sample_metadata_source_gate_20261115.py"
spec=importlib.util.spec_from_file_location("penstemon_sample",P)
m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)

def test_schema_fields_only_existing_source_names():
    sample=[{"name":"sample_accession"},{"name":"scientific_name"},
            {"name":"tax_id"},{"name":"collection_date"}]
    assert m.discover_fields(sample)==["sample_accession","scientific_name","tax_id","collection_date"]

def test_missing_essential_fields_fails_closed():
    with pytest.raises(ValueError,match="Essential"):
        m.discover_fields([{"name":"scientific_name"}])

def test_no_source_taxa_is_metadata_limitation_not_zero_sample():
    raw=("sample_accession\tscientific_name\tsample_alias\nSAMN1\t\t\nSAMN2\t\t\n").encode()
    d=m.count_key_fields(raw)
    assert d["n_unique_samples"]==2
    assert d["original_scientific_name_labels"]=={}
    assert d["original_voucher_join_not_verified"] is True

def test_duplicate_sample_accession_stops():
    raw=b"sample_accession\tscientific_name\nSAMN1\tPenstemon barbatus\nSAMN1\tPenstemon barbatus\n"
    with pytest.raises(ValueError,match="Duplicate"):
        m.count_key_fields(raw)


def test_tabular_returnfields_metadata_parsing_is_not_outcome_access():
    payload="name\tdescription\nsample_accession\tBioSample accession\nscientific_name\tTaxon name\nsample_alias\tOriginal submitter alias\n"
    import json
    # Confirm we can parse source field names without touching biological rows.
    lines=payload.splitlines()
    assert lines[0]=="name\tdescription"
    fields=[{"name":line.split("\t")[0]} for line in lines[1:]]
    assert m.discover_fields(fields)==["sample_accession","sample_alias","scientific_name"]
