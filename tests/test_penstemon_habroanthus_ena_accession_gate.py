"""Synthetic source-only ENA accession / replicate-grain checks."""
import importlib.util
from pathlib import Path
import pytest

P=Path(__file__).resolve().parents[1]/"research"/"prospective"/"penstemon_habroanthus_ena_accession_source_gate_20261115.py"
spec=importlib.util.spec_from_file_location("penstemon_ena",P)
m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)


def test_runs_not_mistaken_for_independent_accessions():
    raw=("run_accession\tsample_accession\tsample_alias\tscientific_name\tstudy_accession\texperiment_accession\n"
       "SRR111\tSAMN1\tp1\tPenstemon barbatus\tSRP1\tSRX1\n"
       "SRR112\tSAMN1\tp1\tPenstemon barbatus\tSRP1\tSRX2\n"
       "SRR113\tSAMN2\tp2\tPenstemon virgatus\tSRP1\tSRX3\n").encode()
    r=m.parse_public_metadata(raw)
    assert r["run_count"]==3
    assert r["sample_accession_count"]==2
    assert r["source_species_label_count"]==2
    assert r["source_species_labels"]["Penstemon barbatus"]==2
    assert r["sample_alias_duplicates"]==["p1"]


def test_empty_public_source_is_not_absence_of_52_manuscript_samples():
    r=m.parse_public_metadata(b"run_accession\tsample_accession\n")
    assert r["status"]=="NO_PUBLIC_RUN_METADATA"
    assert r["not_genome_absence"] is True


def test_duplicate_run_and_missing_sample_fail_closed():
    src=("run_accession\tsample_accession\n"
         "SRR111\tSAMN1\nSRR111\tSAMN2\n").encode()
    with pytest.raises(ValueError,match="Duplicate"):
        m.parse_public_metadata(src)
