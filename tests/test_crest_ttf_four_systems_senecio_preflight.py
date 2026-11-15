"""Offline source-only regression tests; no network and no biological response values."""
import importlib.util
from pathlib import Path

FILE = Path(__file__).resolve().parents[1] / "research" / "prospective" / "ttf_crest_four_systems_senecio_preflight_20261115.py"
spec = importlib.util.spec_from_file_location("senecio_preflight", FILE)
mod = importlib.util.module_from_spec(spec)
spec.loader.exec_module(mod)

def test_source_only_pair_and_pop_keys():
    mock = b"Pair\tEcotype\tPopulation.ID\tFakeOutcome\nD00-H00\tDune\tD00\t999\nD00-H00\tHeadland\tH00\t-999\n"
    actual = mod.parse_source_keys(mock, "phenotype/data_files/allTraitsPairs.txt")
    assert actual["rows"] == 2
    assert actual["population_ids"] == ["D00", "H00"]
    assert actual["pair_counts"] == {"D00-H00": 2}
    assert actual["pair_ecotypes"]["D00-H00"] == ["Dune", "Headland"]

def test_required_keys_fail_closed():
    bad = b"Pair\tEcotype\tFakeOutcome\nD00-H00\tDune\t2\n"
    try:
        mod.parse_source_keys(bad, "phenotype/data_files/allTraitsPairs.txt")
    except ValueError as e:
        assert "metadata columns" in str(e)
    else:
        raise AssertionError("missing population ID must fail")

def test_git_blob_hash_is_sensitive_to_source_changes():
    assert len(mod.git_blob_sha(b"test")) == 40
    assert mod.git_blob_sha(b"test") != mod.git_blob_sha(b"test ")
