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


def test_gbs_imputed_source_provenance_only():
    host = b"Pair,Species1,Species2,GBS_ref,GBS,PctPick1\n1,A,B,0.15897 + nuc * 19.41020,0.7,15\n2,B,C,measured,0.2,99\n"
    mating = b"Pair,nF1M1\n2,10\n"
    result = t.audit_bytes(host, mating, t.sha256_bytes(host), t.sha256_bytes(mating))
    assert result["gbs_provenance_column_present"] is True
    assert result["gbs_proxy_imputed_pair_count"] == 1
    assert result["gbs_provenance_blank_pair_count"] == 0
    assert result["species_labels"] == ["A", "B", "C"]
    assert result["observed_genetic_predictor_admissibility"] == "NOT_QUALIFIED_EXCLUDE_PROXY_GBS"


def test_species_disjoint_geometry_excludes_cross_split_pairs():
    host = b"Pair,Species1,Species2\n1,A,A\n2,B,B\n3,C,C\n4,D,D\n5,A,B\n6,B,C\n"
    mating = b"Pair,nF1M1\n6,1\n1,1\n5,1\n3,1\n4,1\n2,1\n"
    result = t.audit_bytes(host, mating, t.sha256_bytes(host), t.sha256_bytes(mating))
    geom = result["species_disjoint_source_geometry"]
    assert geom["status"] == "GEOMETRY_ONLY_NO_OUTCOME_CLAIM"
    assert geom["n_species"] == 4
    assert geom["n_species_disjoint_nonempty_splits"] > 0
    score = geom["best_balance_capacity"]
    assert set(score["train_species"]).isdisjoint(score["test_species"])
    assert score["train_mating_pair_keys"] > 0
    assert score["test_mating_pair_keys"] > 0
    assert score["cross_species_pair_keys_discarded"] >= 0
