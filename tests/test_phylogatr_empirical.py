from __future__ import annotations

from pathlib import Path

import numpy as np

from ttf.phylogatr_empirical import (
    frozen_edge_mean_p_distances,
    read_aligned_nucleotide_identity,
    sequence_pair_p_distance,
)


def test_pairwise_p_distance_uses_only_jointly_canonical_columns() -> None:
    left = np.frombuffer(b"ACGTNN--", dtype=np.uint8)
    right = np.frombuffer(b"ATGTACGT", dtype=np.uint8)
    assert sequence_pair_p_distance(left, right, minimum_comparable_columns=4) == 0.25
    assert sequence_pair_p_distance(left, right, minimum_comparable_columns=5) is None


def test_frozen_edge_distance_is_mean_over_all_valid_cross_locality_pairs() -> None:
    by_locality = {
        0: (
            np.frombuffer(b"AAAA", dtype=np.uint8),
            np.frombuffer(b"AACC", dtype=np.uint8),
        ),
        1: (
            np.frombuffer(b"AAAT", dtype=np.uint8),
            np.frombuffer(b"CCCC", dtype=np.uint8),
        ),
    }
    result = frozen_edge_mean_p_distances(
        by_locality,
        np.asarray([[0, 1]], dtype=np.int64),
        alignment_length=4,
        minimum_comparable_fraction=0.5,
    )
    # distances: .25, 1.0, .5, .5 -> mean .5625
    assert result.genetic_distance.tolist() == [0.5625]
    assert result.valid_pair_counts.tolist() == [4]
    assert result.minimum_comparable_columns == 2


def test_identity_reader_normalizes_case_but_retains_noncanonical_bytes(tmp_path: Path) -> None:
    path = tmp_path / "aligned.afa"
    path.write_bytes(b">a\nacgtN-\n>b\nACGTry\n")
    alignment = read_aligned_nucleotide_identity(path)
    assert alignment.headers == ("a", "b")
    assert alignment.alignment_length == 6
    assert bytes(alignment.sequences[0]) == b"ACGTN-"
    assert bytes(alignment.sequences[1]) == b"ACGTRY"


def test_identity_reader_preserves_duplicate_header_rows_in_file_order(tmp_path: Path) -> None:
    path = tmp_path / "duplicate.afa"
    path.write_bytes(b">dup\nAAAA\n>dup\nAATT\n>other\nCCCC\n")
    alignment = read_aligned_nucleotide_identity(path)
    assert alignment.headers == ("dup", "dup", "other")
    assert [bytes(x) for x in alignment.sequences] == [b"AAAA", b"AATT", b"CCCC"]
