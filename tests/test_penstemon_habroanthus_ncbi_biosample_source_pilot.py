"""Synthetic NCBI BioSample summary source validation; no real biological data."""
import importlib.util
from pathlib import Path
import json
import pytest

FILE=Path(__file__).resolve().parents[1]/"research"/"prospective"/"penstemon_habroanthus_ncbi_biosample_source_pilot_20261115.py"
spec=importlib.util.spec_from_file_location("penstemon_ncbi_source", FILE)
mod=importlib.util.module_from_spec(spec)
spec.loader.exec_module(mod)


def test_original_summary_yields_only_taxon_and_accession_not_attributes():
    fake={"result":{"uids":["1234","1235"],
      "1234":{"uid":"1234","accession":"SAMN60992931",
          "organism":"Penstemon barbatus","title":"Plant specimen",
          "sampledata":"DO NOT RETURN PRIVATE ORIGINAL ATTRIBUTES"},
      "1235":{"uid":"1235","accession":"SAMN60992932",
          "organism":{"scientificname":"Penstemon virgatus"}}}}
    ret=mod.source_metadata_aggregate(json.dumps(fake).encode(),
        ["SAMN60992931","SAMN60992932"])
    assert ret["status"]=="NCBI_BIOSAMPLE_SOURCE_SUMMARY_ONLY"
    assert ret["result_records"]==2
    assert ret["accessions_with_taxon_identity"]==2
    assert ret["unmatched_source_accessions"]==[]
    assert "DO NOT RETURN" not in repr(ret)
    assert ret["no_sequence_or_genotype_or_phenotype_values_read"] is True


def test_extra_unrequested_accession_fails_closed():
    fake={"result":{"uids":["1234"],
          "1234":{"accession":"SAMN00000001","organism":"Penstemon glaber"}}}
    with pytest.raises(ValueError,match="outside fixed five"):
        mod.source_metadata_aggregate(json.dumps(fake).encode(),["SAMN60992931"])


def test_empty_summary_is_not_negative_biological_result():
    ret=mod.source_metadata_aggregate(
      b'{"result":{"uids":[]}}', ["SAMN60992931"])
    assert ret["result_records"]==0
    assert ret["unmatched_source_accessions"]==["SAMN60992931"]


def test_source_identity_design_is_fixed_before_query():
    assert mod.LIMIT==5
    assert mod.STUDY=="PRJNA1479764"
