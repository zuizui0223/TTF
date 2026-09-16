from __future__ import annotations

import hashlib

import pytest

from ttf.phylogatr_source_projection import project_genes_text_to_animalia


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
