"""Synthetic metadata-only tests; no external Timema outcomes are opened."""
import hashlib
import importlib.util
from pathlib import Path

FILE = Path(__file__).resolve().parents[1] / "research" / "prospective" / "ttf_crest_timema_source_only_20261115.py"
spec = importlib.util.spec_from_file_location("timema_preflight", FILE)
t = importlib.util.module_from_spec(spec)
spec.loader.exec_module(t)


def fixtures():
    host = b"Pair,Species1,Species2,PctPick1,PctPick2,GBS\n1,A,B,0,100,0.9\n2,A,C,0,100,0.2\n3,B,C,0,100,0.4\n"
    mating = b"Pair,nF1M1,yF1M1\n3,123,100\n1,123,10\n"
    return host, mating


def test_source_only_shuffled_pair_keys_and_species_graph():
    host, mating = fixtures()
    d = t.audit_bytes(host, mating, t.sha256_bytes(host), t.sha256_bytes(mating))
    assert d["original_host_pairs"] == 3
    assert d["original_mating_table_pair_keys"] == 2
    assert d["pair_id_overlap"] == 2
    assert d["host_only_pair_ids"] == ["2"]
    assert d["species_graph_component_sizes"] == [3]
    assert d["species_disjoint_transfer_admissibility"] == "NOT_QUALIFIED"


def test_duplicate_pair_keys_rejected():
    host, mating = fixtures()
    duplicate = host + b"2,A,C,0,100,0.2\n"
    try:
        t.audit_bytes(duplicate, mating, t.sha256_bytes(duplicate), t.sha256_bytes(mating))
    except ValueError as e:
        assert "duplicate Pair" in str(e)
    else:
        raise AssertionError("Duplicate Pair was accepted")


def test_orphan_mating_pair_rejected():
    host, _ = fixtures()
    orphan = b"Pair,nF1M1\n4,10\n"
    try:
        t.audit_bytes(host, orphan, t.sha256_bytes(host), t.sha256_bytes(orphan))
    except ValueError as e:
        assert "absent" in str(e)
    else:
        raise AssertionError("Orphan CountData Pair was accepted")


def test_sha_mismatch_is_stop():
    host, mating = fixtures()
    try:
        t.audit_bytes(host, mating, "0" * 64, t.sha256_bytes(mating))
    except ValueError as e:
        assert "SHA256" in str(e)
    else:
        raise AssertionError("source modification was accepted")
