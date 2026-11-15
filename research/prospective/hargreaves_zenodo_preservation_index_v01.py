"""Response-blind Zenodo/Dryad preservation-index preflight for Hargreaves 2020.

Queries official Zenodo public record metadata only. No files or outcome rows
are downloaded. A plausible title or author is NOT enough: source Dryad DOI
must occur as an explicit related/alternate identifier. Search failures or
missing hits are NOT evidence that the archive does not exist.
"""
from __future__ import annotations
import json
import re
from urllib.parse import urlencode
from urllib.request import Request, urlopen

TARGET_DOI = "10.5061/dryad.0vt4b8gv1"
TARGET_FILENAME = "Hargreaves_etal_2020_AmNat_transplant_data.csv"
ENDPOINT = "https://zenodo.org/api/records"
SEARCH_QUERIES = (
    "0vt4b8gv1",
    '"Local adaptation to biotic interactions"',
)
MAX_JSON_BYTES = 750_000
MAX_RECORDS = 25


def doi_norm(value: object) -> str:
    if not isinstance(value, str):
        return ""
    x = value.strip().lower()
    x = re.sub(r"^https?://(?:dx\.)?doi\.org/", "", x)
    x = re.sub(r"^doi:\s*", "", x)
    return x.strip().rstrip("/")


def is_original_dryad_doimatch(rec: dict) -> bool:
    """No title-only, inferred DOI or parent-record ID matches."""
    meta = rec.get("metadata") or {}
    vals = []
    for field in ("identifiers", "related_identifiers", "alternate_identifiers"):
        for item in (meta.get(field) or []):
            if isinstance(item, dict):
                vals.append(item.get("identifier"))
    for field in ("doi", "prereserve_doi"):
        if isinstance(meta.get(field), str):
            vals.append(meta[field])
    for field in ("pids",):
        pid = rec.get(field) or {}
        if isinstance(pid, dict) and isinstance(pid.get("doi"), dict):
            vals.append(pid["doi"].get("identifier"))
    return TARGET_DOI in {doi_norm(v) for v in vals}


def metadata_inventory(rec: dict) -> dict:
    if not is_original_dryad_doimatch(rec):
        raise ValueError("No explicit original Dryad DOI match: no archival identity claim")
    f = rec.get("files") or {}
    if isinstance(f, dict):
        entries = f.get("entries") or {}
        files = list(entries.values()) if isinstance(entries, dict) else []
    elif isinstance(f, list):
        files = f
    else:
        files = []
    names = [str(x.get("key") or x.get("filename") or "") for x in files if isinstance(x, dict)]
    return {
        "record_id": str(rec.get("id", "")),
        "record_url": f"https://zenodo.org/records/{rec.get('id')}",
        "source_original_doi_explicitly_matched": True,
        "file_count_in_public_metadata": len(names),
        "target_original_csv_filename_in_public_metadata": TARGET_FILENAME in names,
        "data_values_or_headers_read": False,
    }


def probe() -> dict:
    attempts, candidates = [], {}
    for query in SEARCH_QUERIES:
        url = ENDPOINT + "?" + urlencode({"q": query, "size": MAX_RECORDS})
        req = Request(url, headers={
            "Accept": "application/json",
            "User-Agent": "TTF-Zenodo-Dryad-source-index-only/2026-11-15",
        })
        try:
            with urlopen(req, timeout=30) as resp:
                raw = resp.read(MAX_JSON_BYTES + 1)
            if len(raw) > MAX_JSON_BYTES:
                raise ValueError("Zenodo search response exceeds source-only metadata cap")
            data = json.loads(raw)
            hits = data.get("hits", {}).get("hits", [])
            if not isinstance(hits, list):
                raise ValueError("Unexpected Zenodo search schema")
            attempts.append({"query": query, "status": "SEARCH_RETURNED",
                             "hits_inspected": len(hits)})
            for record in hits:
                if isinstance(record, dict) and is_original_dryad_doimatch(record):
                    safe = metadata_inventory(record)
                    candidates[safe["record_id"]] = safe
        except Exception as exc:
            attempts.append({"query": query, "status": "ACCESS_ERROR",
                             "error_type": type(exc).__name__,
                             "http_code": getattr(exc, "code", None)})
    verified = sorted(candidates.values(), key=lambda x: x["record_id"])
    return {
        "schema": "ttf_Hargreaves2020_zenodo_official_dryad_preservation_source_only_v01",
        "original_dryad_doi": TARGET_DOI,
        "original_csv_basename": TARGET_FILENAME,
        "official_zenodo_dryad_community": "https://zenodo.org/communities/dryad",
        "source_index_status": ("EXACT_DOI_PUBLIC_PRESERVATION_RECORD_LOCATED"
                                if verified else "NO_EXACT_DOI_PRESERVATION_RECORD_DISCOVERED_BY_QUERIES"),
        "verified_exact_doi_candidates": verified,
        "attempts": attempts,
        "no_candidate_is_NOT_proof_of_absence": True,
        "source_study_overlap_measured": False,
        "independent_new_studies_certified": 0,
        "response_values_read": False,
        "Hargreaves_fitness_outcome_opened": False,
    }


if __name__ == "__main__":
    print("HARGREAVES_ZENODO_SOURCE_PREFLIGHT=" + json.dumps(probe(), sort_keys=True))
