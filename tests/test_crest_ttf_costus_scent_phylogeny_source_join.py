"""Synthetic Costus OSF metadata-only taxon-key and header checks."""
import importlib.util
from pathlib import Path
import pytest

P=Path(__file__).resolve().parents[1]/"research"/"prospective"/"ttf_crest_costus_scent_phylogeny_source_join_20261115.py"
spec=importlib.util.spec_from_file_location("costus_scent_source_join",P)
m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)


def test_only_species_values_are_returned():
    source=b"Species,quantification,pollinator\nCostus allenii,999,Bird\nCostus_lasius,0,Bee\n"
    result=m.inspect_csv_species_keys(source)
    assert result["species_column"]=="Species"
    assert result["row_count"]==2
    assert result["unique_species_labels"]==["costus allenii","costus lasius"]
    assert result["pollinator_label_columns_PRESENT_NOT_INDEPENDENTLY_VERIFIED"]==["pollinator"]
    assert "999" not in str(result)
    assert "Bird" not in str(result)


def test_ambiguous_species_key_does_not_guess():
    source=b"experiment,weight\nCostus allenii,1.2\n"
    result=m.inspect_csv_species_keys(source)
    assert result["species_column"] is None
    assert result["unique_species_labels"] is None
    assert result["status"]=="SOURCE_HEADERS_ONLY_NO_UNAMBIGUOUS_SPECIES_COLUMN"


def test_header_collision_fails_closed():
    source=b"species,species,amount\nA,A,1\n"
    with pytest.raises(ValueError,match="duplicate"):
        m.inspect_csv_species_keys(source)
