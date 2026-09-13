from __future__ import annotations

from dataclasses import dataclass
import hashlib
from pathlib import Path
from typing import Mapping, Sequence

import numpy as np


_CANONICAL_BYTES = frozenset((65, 67, 71, 84, 97, 99, 103, 116))  # A,C,G,T,a,c,g,t


class CharacterMaskError(ValueError):
    pass


@dataclass(frozen=True)
class CanonicalMaskAlignment:
    headers: tuple[str, ...]
    masks: tuple[np.ndarray, ...]
    alignment_length: int

    def by_header(self) -> dict[str, np.ndarray]:
        return {header: mask for header, mask in zip(self.headers, self.masks)}


@dataclass(frozen=True)
class EdgeMaskSupport:
    valid_edges: np.ndarray
    best_comparable_columns: np.ndarray
    minimum_comparable_columns: int

    @property
    def all_edges_valid(self) -> bool:
        return bool(len(self.valid_edges) > 0 and np.all(self.valid_edges))


def _canonical_mask(raw_sequence_line: bytes) -> np.ndarray:
    raw = raw_sequence_line.strip()
    if not raw:
        return np.empty(0, dtype=bool)
    codes = np.frombuffer(raw, dtype=np.uint8)
    # Convert immediately to the only phase-2 information that is allowed to persist.
    mask = (
        (codes == 65)
        | (codes == 67)
        | (codes == 71)
        | (codes == 84)
        | (codes == 97)
        | (codes == 99)
        | (codes == 103)
        | (codes == 116)
    )
    return np.asarray(mask, dtype=bool)


def read_canonical_mask_alignment(path: Path) -> CanonicalMaskAlignment:
    """Read an aligned FASTA into boolean canonical-valid masks only.

    Nucleotide identities are used only transiently to classify each byte as
    canonical A/C/G/T versus noncanonical. The returned object contains no base
    identity and cannot be used to compute nucleotide divergence.
    """
    headers: list[str] = []
    masks: list[np.ndarray] = []
    seen: set[str] = set()
    current_header: str | None = None
    current_chunks: list[np.ndarray] = []

    def finalize() -> None:
        nonlocal current_header, current_chunks
        if current_header is None:
            return
        mask = (
            np.concatenate(current_chunks).astype(bool, copy=False)
            if current_chunks
            else np.empty(0, dtype=bool)
        )
        headers.append(current_header)
        masks.append(mask)
        current_header = None
        current_chunks = []

    with Path(path).open("rb") as handle:
        for raw in handle:
            if raw.startswith(b">"):
                finalize()
                header = raw[1:].strip().decode("utf-8")
                if not header:
                    raise CharacterMaskError("empty aligned FASTA header")
                if header in seen:
                    raise CharacterMaskError(f"duplicate aligned FASTA header: {header}")
                seen.add(header)
                current_header = header
                continue
            if current_header is None:
                if raw.strip():
                    raise CharacterMaskError("sequence content occurs before first FASTA header")
                continue
            chunk = _canonical_mask(raw)
            if len(chunk):
                current_chunks.append(chunk)
    finalize()

    if not headers:
        raise CharacterMaskError("aligned FASTA contains no records")
    lengths = {int(len(mask)) for mask in masks}
    if len(lengths) != 1 or next(iter(lengths)) <= 0:
        raise CharacterMaskError("aligned sequences have unequal or zero lengths")
    alignment_length = int(next(iter(lengths)))
    return CanonicalMaskAlignment(
        headers=tuple(headers),
        masks=tuple(np.asarray(mask, dtype=bool) for mask in masks),
        alignment_length=alignment_length,
    )


def canonical_mask_sha256(alignment: CanonicalMaskAlignment) -> str:
    digest = hashlib.sha256()
    for header, mask in zip(alignment.headers, alignment.masks):
        digest.update(header.encode("utf-8"))
        digest.update(b"\0")
        digest.update(int(len(mask)).to_bytes(8, "big", signed=False))
        digest.update(np.packbits(np.asarray(mask, dtype=np.uint8)).tobytes())
        digest.update(b"\0")
    return digest.hexdigest()


def coordinate_by_header(
    headers: Sequence[str],
    occurrence_rows: Sequence[Mapping[str, str]],
) -> dict[str, tuple[float, float]]:
    by_phylogatr: dict[str, list[Mapping[str, str]]] = {}
    by_accession: dict[str, list[Mapping[str, str]]] = {}
    for row in occurrence_rows:
        pid = str(row.get("phylogatr_id", ""))
        accession = str(row.get("accession", ""))
        if pid:
            by_phylogatr.setdefault(pid, []).append(row)
        if accession:
            by_accession.setdefault(accession, []).append(row)

    out: dict[str, tuple[float, float]] = {}
    for header in headers:
        exact = by_phylogatr.get(str(header), [])
        chosen: Mapping[str, str] | None = None
        if len(exact) == 1:
            chosen = exact[0]
        elif len(exact) == 0:
            fallback = by_accession.get(str(header), [])
            if str(header) and len(fallback) == 1:
                chosen = fallback[0]
        if chosen is None:
            continue
        try:
            lat = float(chosen["latitude"])
            lon = float(chosen["longitude"])
        except (KeyError, TypeError, ValueError):
            continue
        if not np.isfinite(lat) or not np.isfinite(lon):
            continue
        if not -90.0 <= lat <= 90.0 or not -180.0 <= lon <= 180.0:
            continue
        out[str(header)] = (float(lat), float(lon))
    return out


def edge_mask_support(
    masks_by_locality: Mapping[int, Sequence[np.ndarray]],
    edge_nodes: np.ndarray,
    *,
    alignment_length: int,
    minimum_comparable_fraction: float = 0.50,
) -> EdgeMaskSupport:
    nodes = np.asarray(edge_nodes, dtype=np.int64)
    if nodes.ndim != 2 or nodes.shape[1] != 2 or len(nodes) == 0:
        raise ValueError("edge_nodes must be a non-empty m x 2 array")
    if alignment_length < 1:
        raise ValueError("alignment_length must be positive")
    if not 0.0 < float(minimum_comparable_fraction) <= 1.0:
        raise ValueError("minimum_comparable_fraction must be in (0,1]")
    threshold = int(np.ceil(float(minimum_comparable_fraction) * int(alignment_length)))
    valid = np.zeros(len(nodes), dtype=bool)
    best = np.zeros(len(nodes), dtype=np.int64)

    for index, (left, right) in enumerate(nodes):
        left_masks = tuple(masks_by_locality.get(int(left), ()))
        right_masks = tuple(masks_by_locality.get(int(right), ()))
        maximum = 0
        found = False
        for left_mask in left_masks:
            a = np.asarray(left_mask, dtype=bool)
            if len(a) != int(alignment_length):
                raise ValueError("left mask length drift")
            for right_mask in right_masks:
                b = np.asarray(right_mask, dtype=bool)
                if len(b) != int(alignment_length):
                    raise ValueError("right mask length drift")
                overlap = int(np.count_nonzero(a & b))
                maximum = max(maximum, overlap)
                if overlap >= threshold:
                    found = True
                    break
            if found:
                break
        valid[index] = found
        best[index] = int(maximum)
    return EdgeMaskSupport(
        valid_edges=valid,
        best_comparable_columns=best,
        minimum_comparable_columns=threshold,
    )


__all__ = [
    "CanonicalMaskAlignment",
    "CharacterMaskError",
    "EdgeMaskSupport",
    "canonical_mask_sha256",
    "coordinate_by_header",
    "edge_mask_support",
    "read_canonical_mask_alignment",
]
