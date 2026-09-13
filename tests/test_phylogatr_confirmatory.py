from __future__ import annotations

import csv
from pathlib import Path

import numpy as np

from ttf.phylogatr_confirmatory import (
    GENES_HEADERS,
    OCCURRENCE_HEADERS,
    Phase1PanelCandidate,
    choose_one_panel_per_species,
    coordinates_for_headers,
    fasta_header_sha256,
    fasta_headers_only,
    is_coi_family_locus,
    latlon_to_ecef_km,
    phase1_dataset_digest,
    scan_phylogatr_phase1,
)
from ttf.genetic_geometry import prepare_density_scaled_genetic_geometry


ALIASES = (
    "COI",
    "CO1",
    "COX1",
    "COXI",
    "MT-CO1",
    "CYTOCHROME C OXIDASE SUBUNIT I",
    "CYTOCHROME OXIDASE I",
)


def _write_tsv(path: Path, headers: tuple[str, ...], rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=headers, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def _gene_row(
    species: str,
    raw_gene: str,
    raw_dir: str,
    *,
    class_name: str = "Reptilia",
    order: str = "Squamata",
) -> dict[str, str]:
    genus = species.split()[0]
    return {
        "gene": raw_gene,
        "dir": raw_dir,
        "proportion_retained": "1.0",
        "num_seqs_unaligned": "12",
        "num_seqs_aligned": "12",
        "kingdom": "Animalia",
        "phylum": "Chordata",
        "class": class_name,
        "order": order,
        "family": "Exampleidae",
        "genus": genus,
        "species": species,
        "subspecies": "",
        "different_genbank_species": "",
    }


def _write_species_files(
    root: Path,
    raw_dir: str,
    raw_gene: str,
    *,
    n: int = 12,
    sequence_byte: bytes = b"A",
) -> None:
    directory = root / raw_dir
    directory.mkdir(parents=True, exist_ok=True)
    fasta = directory / f"{raw_gene}.afa"
    chunks: list[bytes] = []
    occurrences: list[dict[str, str]] = []
    for index in range(n):
        header = f"ID{index}"
        chunks.extend([f">{header}\n".encode(), sequence_byte * 20 + b"\n"])
        occurrences.append(
            {
                "phylogatr_id": header,
                "accession": f"ACC{index}",
                "source_id": str(index),
                "latitude": str(10.0 + 0.4 * index),
                "longitude": str(20.0 + 0.7 * index),
                "basis_of_record": "PRESERVED_SPECIMEN",
                "coordinate_uncertainty_in_meters": "",
                "issue": "",
                "flag": "",
            }
        )
    fasta.write_bytes(b"".join(chunks))
    _write_tsv(directory / "occurrences.txt", OCCURRENCE_HEADERS, occurrences)


def test_species_prefix_locus_normalization_is_exact() -> None:
    assert is_coi_family_locus("Acanthocercus-annectens-COI", "Acanthocercus annectens", ALIASES)
    assert is_coi_family_locus("Acanthocercus_annectens_MT-CO1", "Acanthocercus annectens", ALIASES)
    assert is_coi_family_locus("COX1", "Acanthocercus annectens", ALIASES)
    assert not is_coi_family_locus("Acanthocercus-annectens-COII", "Acanthocercus annectens", ALIASES)
    assert not is_coi_family_locus("Acanthocercus-annectens-COI-fragment", "Acanthocercus annectens", ALIASES)


def test_fasta_header_reader_never_decodes_nonheader_bytes(tmp_path: Path) -> None:
    path = tmp_path / "x.afa"
    path.write_bytes(b">ID1\n\xff\xfe\xfd\n>ID2\nACGTNN--\n")
    first = fasta_headers_only(path)
    first_hash = fasta_header_sha256(first)
    assert first == ("ID1", "ID2")

    path.write_bytes(b">ID1\nTHIS_SEQUENCE_CHANGED_COMPLETELY\n>ID2\n\x80\x81\x82\n")
    second = fasta_headers_only(path)
    assert second == first
    assert fasta_header_sha256(second) == first_hash


def test_occurrence_join_prefers_unique_phylogatr_id_and_unique_accession() -> None:
    rows = [
        {"phylogatr_id": "H1", "accession": "A1", "latitude": "1", "longitude": "2"},
        {"phylogatr_id": "OTHER", "accession": "H2", "latitude": "3", "longitude": "4"},
        {"phylogatr_id": "D1", "accession": "DUP", "latitude": "5", "longitude": "6"},
        {"phylogatr_id": "D2", "accession": "DUP", "latitude": "7", "longitude": "8"},
    ]
    coordinates = coordinates_for_headers(("H1", "H2", "DUP"), rows)
    np.testing.assert_allclose(coordinates, np.asarray([[1.0, 2.0], [3.0, 4.0]]))


def test_phase1_scan_excludes_aves_chiroptera_decker_and_non_coi(tmp_path: Path) -> None:
    root = tmp_path / "phylogatr-results"
    root.mkdir()
    (root / "cite.txt").write_text("fresh provenance\n")
    rows = [
        _gene_row("Freshus alpha", "Freshus-alpha-COI", "Reptilia/Squamata/Freshus-alpha"),
        _gene_row(
            "Birdus alpha",
            "Birdus-alpha-COI",
            "Aves/Passeriformes/Birdus-alpha",
            class_name="Aves",
            order="Passeriformes",
        ),
        _gene_row(
            "Battus alpha",
            "Battus-alpha-COI",
            "Mammalia/Chiroptera/Battus-alpha",
            class_name="Mammalia",
            order="Chiroptera",
        ),
        _gene_row("Acanthiza nana", "Acanthiza-nana-COI", "Reptilia/Squamata/Acanthiza-nana"),
        _gene_row("Freshus beta", "Freshus-beta-CYTB", "Reptilia/Squamata/Freshus-beta"),
    ]
    _write_tsv(root / "genes.txt", GENES_HEADERS, rows)
    for row in rows:
        _write_species_files(root, row["dir"], row["gene"])

    scan = scan_phylogatr_phase1(
        root,
        aliases=ALIASES,
        excluded_species={"Acanthiza nana"},
        min_localities=12,
        min_endpoint_training_edges=5,
        neighbor_fraction=0.15,
    )
    assert [candidate.species for candidate in scan.candidates] == ["Freshus alpha"]
    statuses = {row["species"]: row["status"] for row in scan.ledger}
    assert statuses["Birdus alpha"] == "excluded_aves"
    assert statuses["Battus alpha"] == "excluded_chiroptera"
    assert statuses["Acanthiza nana"] == "excluded_decker_species"
    assert statuses["Freshus beta"] == "non_coi_locus"


def test_phase1_dataset_digest_is_invariant_to_sequence_identity(tmp_path: Path) -> None:
    root = tmp_path / "phylogatr-results"
    root.mkdir()
    (root / "cite.txt").write_text("fresh provenance\n")
    row = _gene_row("Freshus alpha", "Freshus-alpha-COI", "Reptilia/Squamata/Freshus-alpha")
    _write_tsv(root / "genes.txt", GENES_HEADERS, [row])
    _write_species_files(root, row["dir"], row["gene"], sequence_byte=b"A")

    first_scan = scan_phylogatr_phase1(
        root,
        aliases=ALIASES,
        excluded_species=set(),
        min_localities=12,
        min_endpoint_training_edges=5,
        neighbor_fraction=0.15,
    )
    first_selected = choose_one_panel_per_species(first_scan.candidates)
    first_digest, _ = phase1_dataset_digest(root, first_selected)

    fasta = root / row["dir"] / f"{row['gene']}.afa"
    original_headers = fasta_headers_only(fasta)
    _write_species_files(root, row["dir"], row["gene"], sequence_byte=b"T")
    assert fasta_headers_only(fasta) == original_headers

    second_scan = scan_phylogatr_phase1(
        root,
        aliases=ALIASES,
        excluded_species=set(),
        min_localities=12,
        min_endpoint_training_edges=5,
        neighbor_fraction=0.15,
    )
    second_selected = choose_one_panel_per_species(second_scan.candidates)
    second_digest, _ = phase1_dataset_digest(root, second_selected)
    assert second_digest == first_digest


def _candidate(species: str, n: int, headers: int, gene: str) -> Phase1PanelCandidate:
    latlon = np.column_stack((np.linspace(0.0, 5.0, n), np.linspace(20.0, 30.0, n)))
    ecef = latlon_to_ecef_km(latlon)
    geometry = prepare_density_scaled_genetic_geometry(ecef, neighbor_fraction=0.15)
    # np.unique over ECEF changes row order; the exact order is irrelevant to this tie-break test.
    return Phase1PanelCandidate(
        species=species,
        kingdom="Animalia",
        phylum="Chordata",
        class_name="Reptilia",
        order="Squamata",
        family="Exampleidae",
        genus=species.split()[0],
        raw_gene=gene,
        raw_dir=species.replace(" ", "-"),
        fasta_path=Path(f"{gene}.afa"),
        occurrence_path=Path("occurrences.txt"),
        fasta_headers=tuple(f"H{i}" for i in range(headers)),
        matched_header_count=headers,
        exact_unique_latlon=latlon,
        canonical_latlon=latlon,
        geometry=geometry,
        fasta_header_sha256="x",
        occurrence_sha256="y",
    )


def test_one_panel_rule_prefers_localities_then_header_count_then_gene() -> None:
    panels = [
        _candidate("Freshus alpha", 12, 30, "Freshus-alpha-COX1"),
        _candidate("Freshus alpha", 13, 10, "Freshus-alpha-COI"),
        _candidate("Freshus beta", 12, 10, "Freshus-beta-COX1"),
        _candidate("Freshus beta", 12, 11, "Freshus-beta-MT-CO1"),
    ]
    selected = {panel.species: panel for panel in choose_one_panel_per_species(panels)}
    assert selected["Freshus alpha"].raw_gene == "Freshus-alpha-COI"
    assert selected["Freshus beta"].raw_gene == "Freshus-beta-MT-CO1"
