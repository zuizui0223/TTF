from __future__ import annotations

from collections import Counter, defaultdict
from dataclasses import dataclass
import csv
import hashlib
import json
from pathlib import Path
import re
from typing import Iterable, Mapping, Sequence

import numpy as np

from .genetic_geometry import GeneticSamplingGeometry, prepare_density_scaled_genetic_geometry


EARTH_RADIUS_KM = 6371.0088
GENES_HEADERS = (
    "gene",
    "dir",
    "proportion_retained",
    "num_seqs_unaligned",
    "num_seqs_aligned",
    "kingdom",
    "phylum",
    "class",
    "order",
    "family",
    "genus",
    "species",
    "subspecies",
    "different_genbank_species",
)
OCCURRENCE_HEADERS = (
    "phylogatr_id",
    "accession",
    "source_id",
    "latitude",
    "longitude",
    "basis_of_record",
    "coordinate_uncertainty_in_meters",
    "issue",
    "flag",
)


@dataclass(frozen=True)
class Phase1PanelCandidate:
    species: str
    kingdom: str
    phylum: str
    class_name: str
    order: str
    family: str
    genus: str
    raw_gene: str
    raw_dir: str
    fasta_path: Path
    occurrence_path: Path
    fasta_headers: tuple[str, ...]
    matched_header_count: int
    exact_unique_latlon: np.ndarray
    canonical_latlon: np.ndarray
    geometry: GeneticSamplingGeometry
    fasta_header_sha256: str
    occurrence_sha256: str

    @property
    def n_localities(self) -> int:
        return int(self.geometry.n_localities)

    @property
    def n_headers(self) -> int:
        return int(len(self.fasta_headers))

    def digest_entry(self) -> dict:
        latlon = sorted(
            (float(lat), float(lon)) for lat, lon in np.asarray(self.exact_unique_latlon)
        )
        return {
            "species": self.species,
            "taxonomy": {
                "kingdom": self.kingdom,
                "phylum": self.phylum,
                "class": self.class_name,
                "order": self.order,
                "family": self.family,
                "genus": self.genus,
            },
            "raw_gene": self.raw_gene,
            "raw_dir": self.raw_dir,
            "aligned_header_sha256": self.fasta_header_sha256,
            "occurrence_sha256": self.occurrence_sha256,
            "usable_latlon": [[lat, lon] for lat, lon in latlon],
        }


@dataclass(frozen=True)
class Phase1ScanResult:
    candidates: tuple[Phase1PanelCandidate, ...]
    ledger: tuple[dict, ...]
    genes_rows: int

    @property
    def status_counts(self) -> dict[str, int]:
        return dict(sorted(Counter(str(row["status"]) for row in self.ledger).items()))


def sha256_path(path: Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def sha256_text(value: str) -> str:
    return hashlib.sha256(value.encode("utf-8")).hexdigest()


def canonical_json_sha256(payload: object) -> str:
    encoded = json.dumps(
        payload,
        ensure_ascii=False,
        sort_keys=True,
        separators=(",", ":"),
    ).encode("utf-8")
    return hashlib.sha256(encoded).hexdigest()


def collapse_whitespace(value: str) -> str:
    return " ".join(str(value).strip().split())


def normalize_locus_token(value: str) -> str:
    return collapse_whitespace(re.sub(r"[_\-\s]+", " ", str(value).upper().strip()))


def normalized_locus_after_species_prefix(raw_gene: str, canonical_species: str) -> str:
    gene = normalize_locus_token(raw_gene)
    species = normalize_locus_token(canonical_species)
    prefix = species + " "
    if gene.startswith(prefix):
        gene = gene[len(prefix) :]
    return gene


def is_coi_family_locus(
    raw_gene: str,
    canonical_species: str,
    aliases: Sequence[str],
) -> bool:
    locus = normalized_locus_after_species_prefix(raw_gene, canonical_species)
    normalized_aliases = {normalize_locus_token(alias) for alias in aliases}
    return locus in normalized_aliases


def fasta_headers_only(path: Path) -> tuple[str, ...]:
    """Return FASTA headers while never decoding or retaining sequence lines.

    The file is read in binary mode. Only records whose first byte is ``>`` are
    decoded. Non-header bytes are discarded without character decoding, sequence
    validation, hashing, counting, or any other interpretation.
    """
    headers: list[str] = []
    with Path(path).open("rb") as handle:
        for raw in handle:
            if not raw.startswith(b">"):
                continue
            header = raw[1:].strip().decode("utf-8")
            if header:
                headers.append(header)
    return tuple(headers)


def fasta_header_sha256(headers: Sequence[str]) -> str:
    return sha256_text("\n".join(map(str, headers)) + "\n")


def _read_tsv_exact(path: Path, expected_headers: Sequence[str]) -> list[dict[str, str]]:
    with Path(path).open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if tuple(reader.fieldnames or ()) != tuple(expected_headers):
            raise ValueError(
                f"unexpected TSV schema for {path}: {reader.fieldnames!r} != {tuple(expected_headers)!r}"
            )
        return [{key: str(value or "") for key, value in row.items()} for row in reader]


def read_genes_rows(path: Path) -> list[dict[str, str]]:
    return _read_tsv_exact(path, GENES_HEADERS)


def read_occurrence_rows(path: Path) -> list[dict[str, str]]:
    return _read_tsv_exact(path, OCCURRENCE_HEADERS)


def _valid_coordinate(row: Mapping[str, str]) -> tuple[float, float] | None:
    try:
        lat = float(row["latitude"])
        lon = float(row["longitude"])
    except (KeyError, TypeError, ValueError):
        return None
    if not np.isfinite(lat) or not np.isfinite(lon):
        return None
    if not -90.0 <= lat <= 90.0 or not -180.0 <= lon <= 180.0:
        return None
    return float(lat), float(lon)


def coordinates_for_headers(
    headers: Sequence[str],
    occurrence_rows: Sequence[Mapping[str, str]],
) -> np.ndarray:
    """Join aligned headers to occurrence rows under the frozen exact-match rule."""
    by_phylogatr: dict[str, list[Mapping[str, str]]] = defaultdict(list)
    by_accession: dict[str, list[Mapping[str, str]]] = defaultdict(list)
    for row in occurrence_rows:
        pid = str(row.get("phylogatr_id", ""))
        accession = str(row.get("accession", ""))
        if pid:
            by_phylogatr[pid].append(row)
        if accession:
            by_accession[accession].append(row)

    coordinates: list[tuple[float, float]] = []
    for header in headers:
        exact = by_phylogatr.get(str(header), [])
        chosen: Mapping[str, str] | None = None
        if len(exact) == 1:
            chosen = exact[0]
        elif len(exact) == 0:
            fallback = by_accession.get(str(header), [])
            if str(header) and len(fallback) == 1:
                chosen = fallback[0]
        # Duplicate exact IDs or duplicate accessions are deliberately ambiguous.
        if chosen is None:
            continue
        coordinate = _valid_coordinate(chosen)
        if coordinate is not None:
            coordinates.append(coordinate)
    if not coordinates:
        return np.empty((0, 2), dtype=float)
    return np.asarray(coordinates, dtype=float)


def latlon_to_ecef_km(latlon: np.ndarray) -> np.ndarray:
    x = np.asarray(latlon, dtype=float)
    if x.ndim != 2 or x.shape[1] != 2:
        raise ValueError("latlon must be n x 2")
    lat = np.deg2rad(x[:, 0])
    lon = np.deg2rad(x[:, 1])
    coslat = np.cos(lat)
    return EARTH_RADIUS_KM * np.column_stack(
        [coslat * np.cos(lon), coslat * np.sin(lon), np.sin(lat)]
    )


def canonical_latlon_for_geometry(
    exact_unique_latlon: np.ndarray,
    geometry: GeneticSamplingGeometry,
) -> np.ndarray:
    raw_ecef = latlon_to_ecef_km(exact_unique_latlon)
    lookup = {tuple(map(float, row)): index for index, row in enumerate(raw_ecef)}
    if len(lookup) != len(raw_ecef):
        raise RuntimeError("ECEF conversion collapsed distinct exact lat/lon localities")
    ordered: list[np.ndarray] = []
    for row in np.asarray(geometry.coordinates, dtype=float):
        key = tuple(map(float, row))
        if key not in lookup:
            raise RuntimeError("canonical geometry coordinate missing from source locality map")
        ordered.append(np.asarray(exact_unique_latlon[lookup[key]], dtype=float))
    return np.asarray(ordered, dtype=float)


def _safe_relative_dir(root: Path, raw_dir: str) -> Path:
    rel = Path(str(raw_dir))
    if rel.is_absolute() or ".." in rel.parts:
        raise ValueError(f"unsafe phylogatR dir field: {raw_dir!r}")
    resolved = (Path(root) / rel).resolve()
    root_resolved = Path(root).resolve()
    if root_resolved != resolved and root_resolved not in resolved.parents:
        raise ValueError(f"phylogatR dir escapes root: {raw_dir!r}")
    return resolved


def _candidate_from_row(
    root: Path,
    row: Mapping[str, str],
    *,
    aliases: Sequence[str],
    excluded_species: set[str],
    min_localities: int,
    min_endpoint_training_edges: int,
    neighbor_fraction: float,
) -> tuple[Phase1PanelCandidate | None, str]:
    species = collapse_whitespace(row.get("species", ""))
    if row.get("kingdom", "") != "Animalia":
        return None, "excluded_non_animalia"
    if row.get("class", "") == "Aves":
        return None, "excluded_aves"
    if row.get("order", "") == "Chiroptera":
        return None, "excluded_chiroptera"
    if not species:
        return None, "invalid_species_name"
    if species in excluded_species:
        return None, "excluded_decker_species"
    raw_gene = str(row.get("gene", ""))
    if not is_coi_family_locus(raw_gene, species, aliases):
        return None, "non_coi_locus"

    directory = _safe_relative_dir(root, str(row.get("dir", "")))
    fasta_path = directory / f"{raw_gene}.afa"
    occurrence_path = directory / "occurrences.txt"
    if not fasta_path.is_file():
        return None, "unresolved_aligned_fasta"
    if not occurrence_path.is_file():
        return None, "missing_occurrences"

    headers = fasta_headers_only(fasta_path)
    if not headers:
        return None, "no_aligned_headers"
    occurrence_rows = read_occurrence_rows(occurrence_path)
    latlon = coordinates_for_headers(headers, occurrence_rows)
    if len(latlon) == 0:
        return None, "no_header_occurrence_matches"
    exact_unique = np.unique(latlon, axis=0)
    if len(exact_unique) < int(min_localities):
        return None, "below_min_localities"

    ecef = latlon_to_ecef_km(exact_unique)
    geometry = prepare_density_scaled_genetic_geometry(
        ecef,
        neighbor_fraction=float(neighbor_fraction),
    )
    if geometry.min_endpoint_disjoint_training_edges < int(min_endpoint_training_edges):
        return None, "below_endpoint_disjoint_support"
    canonical_latlon = canonical_latlon_for_geometry(exact_unique, geometry)
    return (
        Phase1PanelCandidate(
            species=species,
            kingdom=str(row.get("kingdom", "")),
            phylum=str(row.get("phylum", "")),
            class_name=str(row.get("class", "")),
            order=str(row.get("order", "")),
            family=str(row.get("family", "")),
            genus=str(row.get("genus", "")),
            raw_gene=raw_gene,
            raw_dir=str(row.get("dir", "")),
            fasta_path=fasta_path,
            occurrence_path=occurrence_path,
            fasta_headers=headers,
            matched_header_count=int(len(latlon)),
            exact_unique_latlon=np.asarray(exact_unique, dtype=float),
            canonical_latlon=np.asarray(canonical_latlon, dtype=float),
            geometry=geometry,
            fasta_header_sha256=fasta_header_sha256(headers),
            occurrence_sha256=sha256_path(occurrence_path),
        ),
        "eligible_candidate",
    )


def scan_phylogatr_phase1(
    root: Path,
    *,
    aliases: Sequence[str],
    excluded_species: Iterable[str],
    min_localities: int = 12,
    min_endpoint_training_edges: int = 5,
    neighbor_fraction: float = 0.15,
) -> Phase1ScanResult:
    root = Path(root)
    genes_path = root / "genes.txt"
    if not genes_path.is_file():
        raise FileNotFoundError(genes_path)
    if not (root / "cite.txt").is_file():
        raise FileNotFoundError(root / "cite.txt")
    rows = read_genes_rows(genes_path)
    exclusion = {collapse_whitespace(name) for name in excluded_species}
    candidates: list[Phase1PanelCandidate] = []
    ledger: list[dict] = []
    for index, row in enumerate(rows):
        candidate, status = _candidate_from_row(
            root,
            row,
            aliases=aliases,
            excluded_species=exclusion,
            min_localities=int(min_localities),
            min_endpoint_training_edges=int(min_endpoint_training_edges),
            neighbor_fraction=float(neighbor_fraction),
        )
        ledger.append(
            {
                "row_index": int(index),
                "species": collapse_whitespace(row.get("species", "")),
                "gene": str(row.get("gene", "")),
                "class": str(row.get("class", "")),
                "order": str(row.get("order", "")),
                "status": status,
            }
        )
        if candidate is not None:
            candidates.append(candidate)
    return Phase1ScanResult(
        candidates=tuple(candidates),
        ledger=tuple(ledger),
        genes_rows=int(len(rows)),
    )


def choose_one_panel_per_species(
    candidates: Sequence[Phase1PanelCandidate],
) -> tuple[Phase1PanelCandidate, ...]:
    grouped: dict[str, list[Phase1PanelCandidate]] = defaultdict(list)
    for candidate in candidates:
        grouped[candidate.species].append(candidate)
    selected: list[Phase1PanelCandidate] = []
    for species in sorted(grouped):
        rows = sorted(
            grouped[species],
            key=lambda row: (
                -row.n_localities,
                -row.n_headers,
                row.raw_gene,
            ),
        )
        selected.append(rows[0])
    return tuple(selected)


def phase1_dataset_digest(
    root: Path,
    selected_before_cap: Sequence[Phase1PanelCandidate],
) -> tuple[str, dict]:
    root = Path(root)
    payload = {
        "cite_sha256": sha256_path(root / "cite.txt"),
        "genes_sha256": sha256_path(root / "genes.txt"),
        "selected_panel_entries": [
            panel.digest_entry() for panel in sorted(selected_before_cap, key=lambda p: p.species)
        ],
    }
    return canonical_json_sha256(payload), payload


def apply_species_cap(
    panels: Sequence[Phase1PanelCandidate],
    *,
    dataset_digest: str,
    maximum_species: int = 250,
) -> tuple[Phase1PanelCandidate, ...]:
    rows = list(panels)
    if len(rows) <= int(maximum_species):
        return tuple(sorted(rows, key=lambda p: p.species))
    keyed = sorted(
        rows,
        key=lambda p: (sha256_text(f"{dataset_digest}|{p.species}"), p.species),
    )
    return tuple(sorted(keyed[: int(maximum_species)], key=lambda p: p.species))


def deterministic_species_split(
    panels: Sequence[Phase1PanelCandidate],
    *,
    dataset_digest: str,
) -> tuple[tuple[str, ...], tuple[str, ...]]:
    keyed = sorted(
        (sha256_text(f"{dataset_digest}|split|{panel.species}"), panel.species)
        for panel in panels
    )
    n_train = len(keyed) // 2
    train = tuple(sorted(species for _, species in keyed[:n_train]))
    evaluation = tuple(sorted(species for _, species in keyed[n_train:]))
    return train, evaluation


__all__ = [
    "EARTH_RADIUS_KM",
    "GENES_HEADERS",
    "OCCURRENCE_HEADERS",
    "Phase1PanelCandidate",
    "Phase1ScanResult",
    "apply_species_cap",
    "canonical_json_sha256",
    "canonical_latlon_for_geometry",
    "choose_one_panel_per_species",
    "collapse_whitespace",
    "coordinates_for_headers",
    "deterministic_species_split",
    "fasta_header_sha256",
    "fasta_headers_only",
    "is_coi_family_locus",
    "latlon_to_ecef_km",
    "normalize_locus_token",
    "normalized_locus_after_species_prefix",
    "phase1_dataset_digest",
    "read_genes_rows",
    "read_occurrence_rows",
    "scan_phylogatr_phase1",
    "sha256_path",
    "sha256_text",
]
