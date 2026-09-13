from __future__ import annotations

from dataclasses import dataclass
import math
from pathlib import Path
from typing import Mapping, Sequence

import numpy as np

from .genetic_geometry import GeneticSamplingGeometry
from .phylogatr_character_mask import coordinate_by_header
from .phylogatr_confirmatory import read_occurrence_rows


_CANONICAL_UPPER = np.asarray([65, 67, 71, 84], dtype=np.uint8)  # A,C,G,T


class EmpiricalSequenceError(ValueError):
    pass


@dataclass(frozen=True)
class AlignedNucleotideIdentity:
    headers: tuple[str, ...]
    sequences: tuple[np.ndarray, ...]
    alignment_length: int

    def by_header(self) -> dict[str, np.ndarray]:
        return {header: sequence for header, sequence in zip(self.headers, self.sequences)}


@dataclass(frozen=True)
class FrozenEdgeGeneticDistances:
    genetic_distance: np.ndarray
    valid_pair_counts: np.ndarray
    minimum_comparable_columns: int
    alignment_length: int


def _upper_ascii(raw: bytes) -> np.ndarray:
    data = np.frombuffer(raw, dtype=np.uint8).copy()
    lower = (data >= 97) & (data <= 122)
    data[lower] -= 32
    return data


def read_aligned_nucleotide_identity(path: Path) -> AlignedNucleotideIdentity:
    """Read aligned FASTA identity after a valid Phase-4 opening authorization.

    The caller is responsible for validating the opening authorization and source
    header hashes *before* invoking this function. Identity is retained in memory
    only; no sequence strings are serialized by this module.
    """
    headers: list[str] = []
    sequences: list[np.ndarray] = []
    seen: set[str] = set()
    current_header: str | None = None
    chunks: list[np.ndarray] = []

    def finalize() -> None:
        nonlocal current_header, chunks
        if current_header is None:
            return
        sequence = (
            np.concatenate(chunks).astype(np.uint8, copy=False)
            if chunks
            else np.empty(0, dtype=np.uint8)
        )
        headers.append(current_header)
        sequences.append(sequence)
        current_header = None
        chunks = []

    with Path(path).open("rb") as handle:
        for raw in handle:
            if raw.startswith(b">"):
                finalize()
                header = raw[1:].strip().decode("utf-8")
                if not header:
                    raise EmpiricalSequenceError("empty aligned FASTA header")
                if header in seen:
                    raise EmpiricalSequenceError(f"duplicate aligned FASTA header: {header}")
                seen.add(header)
                current_header = header
                continue
            if current_header is None:
                if raw.strip():
                    raise EmpiricalSequenceError("sequence content occurs before first FASTA header")
                continue
            chunk = raw.strip()
            if chunk:
                chunks.append(_upper_ascii(chunk))
    finalize()

    if not headers:
        raise EmpiricalSequenceError("aligned FASTA contains no records")
    lengths = {int(len(sequence)) for sequence in sequences}
    if len(lengths) != 1 or next(iter(lengths)) <= 0:
        raise EmpiricalSequenceError("aligned sequences have unequal or zero lengths")
    alignment_length = int(next(iter(lengths)))
    return AlignedNucleotideIdentity(
        headers=tuple(headers),
        sequences=tuple(np.asarray(sequence, dtype=np.uint8) for sequence in sequences),
        alignment_length=alignment_length,
    )


def _canonical_mask(sequence: np.ndarray) -> np.ndarray:
    x = np.asarray(sequence, dtype=np.uint8)
    return (x == 65) | (x == 67) | (x == 71) | (x == 84)


def sequence_pair_p_distance(
    left: np.ndarray,
    right: np.ndarray,
    *,
    minimum_comparable_columns: int,
) -> float | None:
    a = np.asarray(left, dtype=np.uint8)
    b = np.asarray(right, dtype=np.uint8)
    if a.ndim != 1 or b.ndim != 1 or len(a) != len(b):
        raise ValueError("aligned sequences must be equal-length vectors")
    if minimum_comparable_columns < 1:
        raise ValueError("minimum_comparable_columns must be positive")
    comparable = _canonical_mask(a) & _canonical_mask(b)
    n = int(np.count_nonzero(comparable))
    if n < int(minimum_comparable_columns):
        return None
    mismatches = int(np.count_nonzero(a[comparable] != b[comparable]))
    return float(mismatches / n)


def sequences_by_frozen_locality(
    alignment: AlignedNucleotideIdentity,
    occurrence_rows: Sequence[Mapping[str, str]],
    frozen_latlon: np.ndarray,
) -> dict[int, tuple[np.ndarray, ...]]:
    latlon = np.asarray(frozen_latlon, dtype=float)
    if latlon.ndim != 2 or latlon.shape[1] != 2 or len(latlon) == 0:
        raise ValueError("frozen_latlon must be non-empty n x 2")
    lookup = {tuple(map(float, row)): index for index, row in enumerate(latlon)}
    if len(lookup) != len(latlon):
        raise ValueError("frozen_latlon contains duplicate localities")
    coordinates = coordinate_by_header(alignment.headers, occurrence_rows)
    grouped: dict[int, list[np.ndarray]] = {}
    for header, sequence in zip(alignment.headers, alignment.sequences):
        coordinate = coordinates.get(header)
        if coordinate is None:
            continue
        if coordinate not in lookup:
            raise RuntimeError(
                f"authorized sequence coordinate {coordinate!r} is absent from frozen Phase-2 localities"
            )
        index = int(lookup[coordinate])
        grouped.setdefault(index, []).append(np.asarray(sequence, dtype=np.uint8))
    return {index: tuple(values) for index, values in grouped.items()}


def frozen_edge_mean_p_distances(
    sequences_by_locality: Mapping[int, Sequence[np.ndarray]],
    edge_nodes: np.ndarray,
    *,
    alignment_length: int,
    minimum_comparable_fraction: float = 0.50,
) -> FrozenEdgeGeneticDistances:
    nodes = np.asarray(edge_nodes, dtype=np.int64)
    if nodes.ndim != 2 or nodes.shape[1] != 2 or len(nodes) == 0:
        raise ValueError("edge_nodes must be a non-empty m x 2 array")
    if alignment_length < 1:
        raise ValueError("alignment_length must be positive")
    if not 0.0 < float(minimum_comparable_fraction) <= 1.0:
        raise ValueError("minimum_comparable_fraction must lie in (0,1]")
    threshold = int(math.ceil(float(minimum_comparable_fraction) * int(alignment_length)))
    distances = np.empty(len(nodes), dtype=float)
    pair_counts = np.zeros(len(nodes), dtype=np.int64)

    for edge_index, (left, right) in enumerate(nodes):
        left_sequences = tuple(sequences_by_locality.get(int(left), ()))
        right_sequences = tuple(sequences_by_locality.get(int(right), ()))
        valid: list[float] = []
        for a in left_sequences:
            for b in right_sequences:
                distance = sequence_pair_p_distance(
                    a,
                    b,
                    minimum_comparable_columns=threshold,
                )
                if distance is not None:
                    valid.append(distance)
        if not valid:
            raise RuntimeError(
                f"frozen edge {edge_index} unexpectedly lacks a valid cross-locality sequence pair"
            )
        distances[edge_index] = float(np.mean(np.asarray(valid, dtype=float)))
        pair_counts[edge_index] = int(len(valid))

    if not np.isfinite(distances).all() or np.any(distances < 0.0) or np.any(distances > 1.0):
        raise RuntimeError("empirical locality-pair p-distance is outside [0,1]")
    return FrozenEdgeGeneticDistances(
        genetic_distance=distances,
        valid_pair_counts=pair_counts,
        minimum_comparable_columns=threshold,
        alignment_length=int(alignment_length),
    )


def extract_species_frozen_edge_distances(
    fasta_path: Path,
    occurrence_path: Path,
    frozen_latlon: np.ndarray,
    geometry: GeneticSamplingGeometry,
    *,
    minimum_comparable_fraction: float = 0.50,
) -> FrozenEdgeGeneticDistances:
    alignment = read_aligned_nucleotide_identity(fasta_path)
    occurrences = read_occurrence_rows(occurrence_path)
    grouped = sequences_by_frozen_locality(alignment, occurrences, frozen_latlon)
    return frozen_edge_mean_p_distances(
        grouped,
        geometry.edge_nodes,
        alignment_length=alignment.alignment_length,
        minimum_comparable_fraction=minimum_comparable_fraction,
    )


__all__ = [
    "AlignedNucleotideIdentity",
    "EmpiricalSequenceError",
    "FrozenEdgeGeneticDistances",
    "extract_species_frozen_edge_distances",
    "frozen_edge_mean_p_distances",
    "read_aligned_nucleotide_identity",
    "sequence_pair_p_distance",
    "sequences_by_frozen_locality",
]
