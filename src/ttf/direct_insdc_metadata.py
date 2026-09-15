"""Response-blind GenBank metadata extraction for the direct INSDC arm.

This module deliberately never parses ORIGIN sequence content. It consumes GenBank
flatfile records and emits only source/annotation metadata needed to construct a
candidate geometry panel before nucleotide identity may be opened.
"""
from __future__ import annotations

from dataclasses import dataclass
import re
from typing import Iterable, Iterator

_LATLON_RE = re.compile(r'^\s*/lat_lon="?([^"\n]+)"?\s*$')
_ORGANISM_RE = re.compile(r'^\s+ORGANISM\s+(.+?)\s*$')
_ACCESSION_RE = re.compile(r'^ACCESSION\s+(\S+)')
_LOCUS_TERMS = ("COI", "CO1", "COX1", "CYTOCHROME C OXIDASE SUBUNIT I", "CYTOCHROME OXIDASE SUBUNIT I")

@dataclass(frozen=True)
class InsdcMetadataRecord:
    accession: str
    organism: str
    taxonomy: tuple[str, ...]
    lat_lon: str
    coi_family_annotated: bool


def _coi_annotation(text: str) -> bool:
    upper = text.upper()
    return any(term in upper for term in _LOCUS_TERMS)


def iter_genbank_metadata(lines: Iterable[str]) -> Iterator[InsdcMetadataRecord]:
    """Yield metadata-only records without interpreting sequence lines.

    ORIGIN starts a hard ignore state lasting until the record terminator `//`.
    This is a firewall, not an optimization: sequence characters are never added
    to the annotation buffer and therefore cannot influence eligibility.
    """
    accession = ""
    organism = ""
    taxonomy_parts: list[str] = []
    lat_lon = ""
    annotation_lines: list[str] = []
    in_origin = False
    after_organism = False

    def emit():
        if accession and organism and lat_lon:
            text = "\n".join(annotation_lines)
            return InsdcMetadataRecord(
                accession=accession,
                organism=organism,
                taxonomy=tuple(taxonomy_parts),
                lat_lon=lat_lon,
                coi_family_annotated=_coi_annotation(text),
            )
        return None

    for raw in lines:
        line = raw.rstrip("\n")
        if line.startswith("ORIGIN"):
            in_origin = True
            continue
        if line.startswith("//"):
            rec = emit()
            if rec is not None:
                yield rec
            accession = organism = lat_lon = ""
            taxonomy_parts = []
            annotation_lines = []
            in_origin = False
            after_organism = False
            continue
        if in_origin:
            continue
        m = _ACCESSION_RE.match(line)
        if m:
            accession = m.group(1)
        m = _ORGANISM_RE.match(line)
        if m:
            organism = m.group(1).strip()
            after_organism = True
            continue
        if after_organism:
            if line.startswith("REFERENCE") or line.startswith("FEATURES"):
                after_organism = False
            elif line.startswith("            "):
                taxonomy_parts.extend(x.strip() for x in line.strip().rstrip(".").split(";") if x.strip())
        m = _LATLON_RE.match(line)
        if m:
            lat_lon = m.group(1).strip()
        annotation_lines.append(line)
