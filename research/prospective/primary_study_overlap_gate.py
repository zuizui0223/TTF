"""Response-blind bibliographic overlap gate for new TTF transplant sources.

Only accepts a minimal metadata-only CSV. It cannot consume original fitness
CSVs, infer independent replication from distinct dataset DOIs, or certify
that different papers used different experiments.
"""
from __future__ import annotations

import argparse
import csv
import json
import re
from collections import defaultdict
from pathlib import Path

FIELDS = (
    "source_dataset_doi", "source_study_id", "primary_paper_doi",
    "primary_first_author", "primary_year", "primary_title", "plant_taxon",
)
MAX_ROWS = 100000


def norm_doi(value: str) -> str:
    value = value.strip().lower()
    value = re.sub(r"^(?:https?://(?:dx\.)?doi\.org/|doi:\s*)", "", value)
    if not value:
        return ""
    if not re.fullmatch(r"10\.\d{4,9}/\S+", value):
        raise ValueError("Malformed DOI in metadata-only manifest")
    return value.rstrip("./")


def norm_text(value: str) -> str:
    return " ".join(re.findall(r"[a-z0-9]+", value.casefold()))


def read_manifest(path: Path) -> dict:
    """Fail closed on response-bearing or extra columns, even if unused."""
    if path.stat().st_size > 15_000_000:
        raise ValueError("Metadata manifest exceeds bounded input size")
    with path.open("r", encoding="utf-8-sig", newline="") as stream:
        reader = csv.DictReader(stream, restkey="__OVERFLOW__", strict=True)
        if reader.fieldnames is None or tuple(reader.fieldnames) != FIELDS:
            raise ValueError("Only exact metadata-only header is authorized; refuse outcome columns")
        rows = list(reader)
    studies = {}
    source_id = None
    if len(rows) > MAX_ROWS:
        raise ValueError("Metadata manifest exceeds bounded row count")
    for row in rows:
        if "__OVERFLOW__" in row or any(v is None for v in row.values()):
            raise ValueError("Malformed or oversized metadata manifest")
        dataset = norm_doi(row["source_dataset_doi"])
        sid = row["source_study_id"].strip()
        doi = norm_doi(row["primary_paper_doi"])
        author = norm_text(row["primary_first_author"])
        year = row["primary_year"].strip()
        title = norm_text(row["primary_title"])
        taxon = norm_text(row["plant_taxon"])
        if not dataset or not sid or not taxon or not re.fullmatch(r"(?:18|19|20)\d{2}", year):
            raise ValueError("Missing mandatory dataset/study/taxon or invalid year")
        if source_id is None:
            source_id = dataset
        if dataset != source_id:
            raise ValueError("One metadata manifest must represent one source dataset")
        if not doi and (not title or not author):
            raise ValueError("Primary DOI or complete title/author information required")
        ident = (doi, author, year, title)
        if sid in studies:
            if studies[sid]["identity"] != ident:
                raise ValueError("Inconsistent bibliographic identity for original study")
            studies[sid]["taxa"].add(taxon)
        else:
            studies[sid] = {"identity": ident, "taxa": {taxon}}
    if not studies:
        raise ValueError("Empty source-only metadata manifest")
    return {"source_dataset_doi": source_id, "studies": studies, "rows": len(rows)}


def audit_overlap(base: dict, candidate: dict) -> dict:
    if base["source_dataset_doi"] == candidate["source_dataset_doi"]:
        raise ValueError("Cannot treat the same source dataset as independent")
    known = defaultdict(set)
    for key, entry in base["studies"].items():
        doi, author, year, title = entry["identity"]
        if doi:
            known[("doi", doi)].add(key)
        if title:
            known[("title-year", title, year)].add(key)
        if author:
            known[("author-year", author, year)].add(key)
    exact = probable = unmatched = 0
    details = []
    for sid, entry in sorted(candidate["studies"].items()):
        doi, author, year, title = entry["identity"]
        doi_refs = known.get(("doi", doi), set()) if doi else set()
        title_refs = known.get(("title-year", title, year), set()) if title else set()
        author_refs = known.get(("author-year", author, year), set()) if author else set()
        if doi_refs or title_refs:
            status = "EXACT_BIBLIOGRAPHIC_OVERLAP"
            exact += 1
        elif author_refs:
            status = "POSSIBLE_SAME_PRIMARY_STUDY_REVIEW"
            probable += 1
        else:
            status = "UNLINKED_NOT_CERTIFIED_INDEPENDENT"
            unmatched += 1
        details.append({"candidate_source_study_id": sid, "status": status})
    return {
        "schema": "ttf_source_only_primary_study_overlap_audit_v01",
        "reference_dataset_doi": base["source_dataset_doi"],
        "candidate_dataset_doi": candidate["source_dataset_doi"],
        "reference_primary_studies": len(base["studies"]),
        "candidate_primary_studies": len(candidate["studies"]),
        "confirmed_bibliographic_overlap": exact,
        "possible_primary_study_overlap_needs_review": probable,
        "unlinked_but_not_independence_certified": unmatched,
        "verified_independent_primary_studies": 0,
        "outcome_values_read": False,
        "distinct_dataset_DOIs_do_NOT_guarantee_independent_experiments": True,
        "candidate_results": details,
        "status": "SOURCE_ONLY_MANUAL_PRIMARY_EXPERIMENT_ADJUDICATION_REQUIRED",
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--reference", type=Path, required=True)
    parser.add_argument("--candidate", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    output = audit_overlap(read_manifest(args.reference), read_manifest(args.candidate))
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(output, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(output["status"])


if __name__ == "__main__":
    main()
