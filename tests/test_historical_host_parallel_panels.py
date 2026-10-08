import csv
import hashlib

import pytest

from scripts.audit_historical_host_parallel_panels import (
    compute_overlap,
    read_pinned_panel,
)


def _write(path, pairs):
    with path.open("w", newline="", encoding="utf-8") as handle:
        w = csv.DictWriter(handle, fieldnames=["species", "accepted_host_ids"], lineterminator="\n")
        w.writeheader()
        w.writerows([{"species": name, "accepted_host_ids": ids} for name, ids in pairs])
    return hashlib.sha256(path.read_bytes()).hexdigest()


def test_shared_species_exact_same_hosts_independent_of_id_format(tmp_path):
    one = tmp_path / "a.csv"
    two = tmp_path / "b.csv"
    sha1 = _write(one, [("A a", "1.0;2.0"), ("B b", "3.0")])
    sha2 = _write(two, [("A a", "2;1"), ("C c", "4")])
    first = read_pinned_panel(one, sha1, 2)
    second = read_pinned_panel(two, sha2, 2)
    got = compute_overlap(first, second)
    assert got["overlap_species"] == 1
    assert got["union_species"] == 3
    assert got["overlap_accepted_host_ids_identical_species"] == 1
    assert got["overlap_accepted_host_id_disagreements"] == []


def test_overlap_detects_true_diet_source_disagreement(tmp_path):
    p1, p2 = tmp_path / "a.csv", tmp_path / "b.csv"
    sha1 = _write(p1, [("A a", "1;2")])
    sha2 = _write(p2, [("A a", "2;3")])
    got = compute_overlap(read_pinned_panel(p1, sha1, 1), read_pinned_panel(p2, sha2, 1))
    assert got["overlap_accepted_host_id_disagreements"] == ["A a"]


def test_frozen_byte_drift_is_not_rescued_by_matching_species(tmp_path):
    p = tmp_path / "c.csv"
    sha = _write(p, [("A a", "1")])
    p.write_bytes(p.read_bytes().replace(b"1\n", b"1.0\n"))
    with pytest.raises(RuntimeError, match="SHA256 mismatch"):
        read_pinned_panel(p, sha, 1)


def test_duplicate_species_and_host_ids_are_rejected(tmp_path):
    for pairs in ([("A a", "1"), ("A a", "1")], [("A a", "1;1.0")]):
        p = tmp_path / "d.csv"
        sha = _write(p, pairs)
        with pytest.raises(RuntimeError, match="invalid or duplicate"):
            read_pinned_panel(p, sha, len(pairs))
