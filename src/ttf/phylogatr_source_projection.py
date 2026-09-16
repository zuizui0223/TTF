from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
import hashlib
from pathlib import Path

from .phylogatr_confirmatory import GENES_HEADERS


@dataclass(frozen=True)
class GenesProjectionResult:
    projected_text: str
    original_rows: int
    retained_rows: int
    removed_rows: int
    kingdom_counts: dict[str, int]
    original_sha256: str
    projected_sha256: str
    sequence_identity_opened: bool = False


def _sha256_text(value: str) -> str:
    return hashlib.sha256(value.encode("utf-8")).hexdigest()


def _sha256_bytes(value: bytes) -> str:
    return hashlib.sha256(value).hexdigest()


def project_genes_text_to_animalia(text: str) -> GenesProjectionResult:
    """Project phylogatR genes metadata to Animalia without touching sequence data.

    The header and every retained source row are copied byte-for-byte at the
    decoded UTF-8 text level, including their original line endings. Only rows
    whose exact ``kingdom`` field is ``Animalia`` are retained.
    """
    if not isinstance(text, str) or not text:
        raise ValueError("genes.txt must be non-empty UTF-8 text")
    lines = text.splitlines(keepends=True)
    if not lines:
        raise ValueError("genes.txt must contain a header")
    header_fields = tuple(lines[0].rstrip("\r\n").split("\t"))
    if header_fields != tuple(GENES_HEADERS):
        raise ValueError("genes.txt schema does not match the frozen phylogatR schema")

    kingdom_index = GENES_HEADERS.index("kingdom")
    retained = [lines[0]]
    counts: Counter[str] = Counter()
    original_rows = 0
    retained_rows = 0
    for line in lines[1:]:
        fields = line.rstrip("\r\n").split("\t")
        if len(fields) != len(GENES_HEADERS):
            raise ValueError("genes.txt schema drift in data row")
        original_rows += 1
        kingdom = fields[kingdom_index]
        counts[kingdom] += 1
        if kingdom == "Animalia":
            retained.append(line)
            retained_rows += 1

    projected = "".join(retained)
    return GenesProjectionResult(
        projected_text=projected,
        original_rows=original_rows,
        retained_rows=retained_rows,
        removed_rows=original_rows - retained_rows,
        kingdom_counts=dict(sorted(counts.items())),
        original_sha256=_sha256_text(text),
        projected_sha256=_sha256_text(projected),
        sequence_identity_opened=False,
    )


def project_genes_file_in_place(root: Path) -> dict[str, object]:
    """Replace only ``genes.txt`` with its deterministic Animalia projection.

    ``cite.txt`` and every per-species file remain untouched. The returned
    receipt records both pre-projection and projected metadata hashes and makes
    the response-blind firewall explicit.
    """
    root = Path(root)
    genes_path = root / "genes.txt"
    cite_path = root / "cite.txt"
    if not genes_path.is_file() or not cite_path.is_file():
        raise FileNotFoundError("phylogatR root must contain genes.txt and cite.txt")

    raw_genes = genes_path.read_bytes()
    try:
        text = raw_genes.decode("utf-8")
    except UnicodeDecodeError as exc:
        raise ValueError("genes.txt must be valid UTF-8") from exc
    projection = project_genes_text_to_animalia(text)
    projected_bytes = projection.projected_text.encode("utf-8")
    genes_path.write_bytes(projected_bytes)

    return {
        "schema": "ttf_genetic_phylogatr_source_projection_receipt_v0.1",
        "status": "ANIMALIA_METADATA_PROJECTION_APPLIED",
        "raw_genes_sha256": _sha256_bytes(raw_genes),
        "projected_genes_sha256": _sha256_bytes(projected_bytes),
        "cite_sha256": _sha256_bytes(cite_path.read_bytes()),
        "original_row_count": projection.original_rows,
        "retained_animalia_row_count": projection.retained_rows,
        "removed_non_animalia_row_count": projection.removed_rows,
        "kingdom_counts": projection.kingdom_counts,
        "sequence_identity_opened": False,
        "pairwise_genetic_distances_opened": False,
        "empirical_ttf_opened": False,
    }


__all__ = [
    "GenesProjectionResult",
    "project_genes_file_in_place",
    "project_genes_text_to_animalia",
]
