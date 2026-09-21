import importlib.util
import json
from pathlib import Path

SCRIPT = Path("scripts/freeze_relational_fresh_candidate_census.py")


def load():
    spec = importlib.util.spec_from_file_location("relcensus", SCRIPT)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_species_hash_is_deterministic_source_bound():
    m = load()
    assert m.species_hash("a" * 64, "Alpha beta") == m.species_hash("a" * 64, "Alpha beta")
    assert m.species_hash("a" * 64, "Alpha beta") != m.species_hash("b" * 64, "Alpha beta")


def test_read_exclusion_species_phase1_manifest(tmp_path):
    m = load()
    p = tmp_path / "phase1.json"
    p.write_text(json.dumps({"selected_panels": {"Beta gamma": {}, "Alpha beta": {}}}))
    assert m.read_exclusion_species(p) == {"Alpha beta", "Beta gamma"}


def test_read_exclusion_species_csv(tmp_path):
    m = load()
    p = tmp_path / "prior.csv"
    p.write_text("species,value\nAlpha beta,1\nBeta gamma,2\n")
    assert m.read_exclusion_species(p) == {"Alpha beta", "Beta gamma"}


def test_species_digest_preserves_frozen_order():
    m = load()
    assert m.species_digest(["Alpha beta", "Beta gamma"]) != m.species_digest(["Beta gamma", "Alpha beta"])
