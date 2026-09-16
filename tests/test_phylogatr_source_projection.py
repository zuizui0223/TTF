from __future__ import annotations

import hashlib
from pathlib import Path

import pytest

from ttf.phylogatr_source_projection import (
    project_genes_file_in_place,
    project_genes_text_to_animalia,
)


HEADER = (
    "gene\tdir\tproportion_retained\tnum_seqs_unaligned\tnum_seqs_aligned\tkingdom\t"
    "phylum\tclass\torder\tfamily\tgenus\tspecies\tsubspecies\tdifferent_genbank_species\n"
)


def _row(*, gene: str, kingdom: str, species: str) -> str:
    return (
        f"{gene}\tInsecta/Test/Testidae/{species.replace(' ', '-')}\t1.0\t3\t3\t{kingdom}\t"
        f"Arthropoda\tInsecta\tTestOrder\tTestidae\tTestgenus\t{species}\t\t\n"
    )


def test_projection_keeps_only_animalia_rows_byte_stably() -> None:
    animal_a = _row(gene="Alpha-beta-COI", kingdom="Animalia", species="Alpha beta")
    plant = _row(gene="Planta-gamma-COI", kingdom="Plantae", species="Planta gamma")
    animal_b = _row(gene="Delta-epsilon-CYTB", kingdom="Animalia", species="Delta epsilon")
    raw = HEADER + animal_a + plant + animal_b

    result = project_genes_text_to_animalia(raw)

    expected = HEADER + animal_a + animal_b
    assert result.projected_text == expected
    assert result.original_rows == 3
    assert result.retained_rows == 2
    assert result.removed_rows == 1
    assert result.kingdom_counts == {"Animalia": 2, "Plantae": 1}
    assert result.original_sha256 == hashlib.sha256(raw.encode("utf-8")).hexdigest()
    assert result.projected_sha256 == hashlib.sha256(expected.encode("utf-8")).hexdigest()
    assert result.sequence_identity_opened is False


def test_projection_rejects_schema_drift() -> None:
    malformed = "gene\tkingdom\nAlpha-beta-COI\tAnimalia\n"
    with pytest.raises(ValueError, match="genes.txt schema"):
        project_genes_text_to_animalia(malformed)


def test_file_projection_rewrites_only_genes_and_returns_provenance(tmp_path: Path) -> None:
    root = tmp_path / "phylogatr-results"
    root.mkdir()
    animal = _row(gene="Alpha-beta-COI", kingdom="Animalia", species="Alpha beta")
    plant = _row(gene="Planta-gamma-COI", kingdom="Plantae", species="Planta gamma")
    raw = HEADER + animal + plant
    genes = root / "genes.txt"
    cite = root / "cite.txt"
    genes.write_text(raw, encoding="utf-8")
    cite.write_bytes(b"source provenance\n")

    receipt = project_genes_file_in_place(root)

    expected = HEADER + animal
    assert genes.read_text(encoding="utf-8") == expected
    assert cite.read_bytes() == b"source provenance\n"
    assert receipt["schema"] == "ttf_genetic_phylogatr_source_projection_receipt_v0.1"
    assert receipt["raw_genes_sha256"] == hashlib.sha256(raw.encode()).hexdigest()
    assert receipt["projected_genes_sha256"] == hashlib.sha256(expected.encode()).hexdigest()
    assert receipt["cite_sha256"] == hashlib.sha256(b"source provenance\n").hexdigest()
    assert receipt["original_row_count"] == 2
    assert receipt["retained_animalia_row_count"] == 1
    assert receipt["removed_non_animalia_row_count"] == 1
    assert receipt["kingdom_counts"] == {"Animalia": 1, "Plantae": 1}
    assert receipt["sequence_identity_opened"] is False
