from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from ttf.phylogatr_character_mask import (
    CharacterMaskError,
    canonical_mask_sha256,
    edge_mask_support,
    read_canonical_mask_alignment,
)


def _write(path: Path, records: list[tuple[str, str]]) -> None:
    lines: list[str] = []
    for header, sequence in records:
        lines.extend([f">{header}\n", sequence + "\n"])
    path.write_text("".join(lines))


def test_character_identity_is_erased_from_mask_and_hash(tmp_path: Path) -> None:
    a = tmp_path / "a.afa"
    c = tmp_path / "c.afa"
    _write(a, [("H1", "AAAAacgt"), ("H2", "TTTTGCAT")])
    _write(c, [("H1", "CCCCtgca"), ("H2", "GGGGATGC")])
    aa = read_canonical_mask_alignment(a)
    cc = read_canonical_mask_alignment(c)
    assert aa.headers == cc.headers == ("H1", "H2")
    for left, right in zip(aa.masks, cc.masks):
        assert np.array_equal(left, right)
        assert np.all(left)
    assert canonical_mask_sha256(aa) == canonical_mask_sha256(cc)


def test_noncanonical_gap_n_and_iupac_are_invalid_but_lowercase_acgt_valid(tmp_path: Path) -> None:
    path = tmp_path / "mask.afa"
    _write(path, [("H1", "ACGTacgtNRY-?.")])
    alignment = read_canonical_mask_alignment(path)
    expected = np.asarray(
        [True, True, True, True, True, True, True, True, False, False, False, False, False, False],
        dtype=bool,
    )
    assert np.array_equal(alignment.masks[0], expected)


def test_mask_hash_changes_when_valid_invalid_pattern_changes(tmp_path: Path) -> None:
    first = tmp_path / "first.afa"
    second = tmp_path / "second.afa"
    _write(first, [("H", "ACGTNN")])
    _write(second, [("H", "ACG-NN")])
    assert canonical_mask_sha256(read_canonical_mask_alignment(first)) != canonical_mask_sha256(
        read_canonical_mask_alignment(second)
    )


def test_duplicate_header_is_rejected(tmp_path: Path) -> None:
    path = tmp_path / "dup.afa"
    _write(path, [("H", "AAAA"), ("H", "CCCC")])
    with pytest.raises(CharacterMaskError, match="duplicate aligned FASTA header"):
        read_canonical_mask_alignment(path)


def test_unequal_or_zero_alignment_length_is_rejected(tmp_path: Path) -> None:
    unequal = tmp_path / "unequal.afa"
    _write(unequal, [("H1", "AAAA"), ("H2", "AAA")])
    with pytest.raises(CharacterMaskError, match="unequal or zero"):
        read_canonical_mask_alignment(unequal)

    zero = tmp_path / "zero.afa"
    _write(zero, [("H1", ""), ("H2", "")])
    with pytest.raises(CharacterMaskError, match="unequal or zero"):
        read_canonical_mask_alignment(zero)


def test_wrapped_sequence_lines_are_concatenated_before_masking(tmp_path: Path) -> None:
    path = tmp_path / "wrapped.afa"
    path.write_bytes(b">H1\nACGT\nNN--\n>H2\nacgt\nRY??\n")
    alignment = read_canonical_mask_alignment(path)
    assert alignment.alignment_length == 8
    assert np.array_equal(
        alignment.masks[0],
        np.asarray([True, True, True, True, False, False, False, False]),
    )
    assert np.array_equal(alignment.masks[1], alignment.masks[0])


def test_edge_support_uses_exact_fifty_percent_ceiling() -> None:
    # Alignment length 5 -> ceil(2.5)=3 comparable canonical columns.
    edges = np.asarray([[0, 1], [1, 2]], dtype=np.int64)
    masks = {
        0: [np.asarray([1, 1, 1, 0, 0], dtype=bool)],
        1: [np.asarray([1, 1, 1, 1, 0], dtype=bool)],
        2: [np.asarray([1, 1, 0, 0, 0], dtype=bool)],
    }
    support = edge_mask_support(
        masks,
        edges,
        alignment_length=5,
        minimum_comparable_fraction=0.50,
    )
    assert support.minimum_comparable_columns == 3
    assert support.valid_edges.tolist() == [True, False]
    assert support.best_comparable_columns.tolist() == [3, 2]
    assert not support.all_edges_valid


def test_edge_support_accepts_any_one_cross_locality_pair() -> None:
    edges = np.asarray([[0, 1]], dtype=np.int64)
    masks = {
        0: [
            np.asarray([1, 0, 0, 0], dtype=bool),
            np.asarray([1, 1, 1, 0], dtype=bool),
        ],
        1: [
            np.asarray([1, 1, 0, 0], dtype=bool),
            np.asarray([1, 1, 1, 1], dtype=bool),
        ],
    }
    support = edge_mask_support(
        masks,
        edges,
        alignment_length=4,
        minimum_comparable_fraction=0.50,
    )
    assert support.minimum_comparable_columns == 2
    assert support.valid_edges.tolist() == [True]
    assert support.best_comparable_columns.tolist() == [3]
    assert support.all_edges_valid
