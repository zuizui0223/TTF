"""Official Zenodo Dryad preservation metadata search for Bontrager 2021.

Uses public Zenodo record metadata ONLY (no archive members or fitness rows).
A candidate must explicitly reference the exact original Dryad DOI. Title and
author similarity alone never authorize treating files as the source original.
"""
from __future__ import annotations
import json
import re
from urllib.parse import urlencode
from urllib.request import Request, urlopen
from urllib.error import HTTPError, URLError

TARGET_DOI = "10.5061/dryad.pnvx0k6mc"
TARGET_FILENAME = "dryad.zip"
INDEX_URL = "https://zenodo.org/api/records"
QUERIES = ("pnvx0k6mc", '"Adaptation across geographic ranges is consistent with strong selection"')
MAX_BYTES = 750_000
MAX_HITS = 25

def canonical_doi(x: object) -> str:
    if not isinstance(x, str):
        return ""
    return re.sub(r"^(?:https?://(?:dx\.)?doi\.org/|doi:)", "", x.strip().lower()).strip().rstrip("/")

def exact_original_doi_in_record(record: dict) -> bool:
    meta = record.get("metadata") or {}
    identifiers = []
    for f in ("related_identifiers", "identifiers", "alternate_identifiers"):
        for item in (meta.get(f) or []):
            if isinstance(item, dict):
                identifiers.append(item.get("identifier"))
    identifiers.extend([meta.get("doi"), record.get("doi")])
    pids = record.get("pids") or {}
    if isinstance(pids.get("doi"), dict):
        identifiers.append(pids["doi"].get("identifier"))
    return TARGET_DOI in {canonical_doi(value) for value in identifiers}

def original_public_files(record: dict) -> list:
    files = record.get("files") or {}
    if isinstance(files, list):
        data = files
    elif isinstance(files, dict) and isinstance(files.get("entries"), dict):
        data = list(files["entries"].values())
    else:
        data = []
    return [{
        "filename": item.get("key") or item.get("filename"),
        "byte_length": item.get("size"),
        "published_checksum": item.get("checksum"),
    } for item in data if isinstance(item, dict)]

def probe() -> dict:
    attempts = []
    matches = {}
    for query in QUERIES:
        url = INDEX_URL + "?" + urlencode({"q": query, "size": MAX_HITS})
        req = Request(url, headers={"Accept": "application/json",
                      "User-Agent": "TTF-Bontrager-original-Zenodo-index-only/2026-11-15"})
        try:
            with urlopen(req, timeout=35) as response:
                buf = response.read(MAX_BYTES + 1)
            if len(buf) > MAX_BYTES:
                raise ValueError("Public metadata response cap exceeded")
            hits = json.loads(buf).get("hits", {}).get("hits", [])
            if not isinstance(hits, list):
                raise ValueError("Unexpected Zenodo API result type")
            attempts.append({"query": query, "result": "SEARCH_RETURNED", "records": len(hits)})
            for r in hits:
                if not isinstance(r, dict) or not exact_original_doi_in_record(r):
                    continue
                members = original_public_files(r)
                rid = str(r.get("id",""))
                matches[rid] = {
                    "zenodo_record_id": rid,
                    "zenodo_public_url": "https://zenodo.org/records/" + rid,
                    "original_dryad_doi_explicit_exact_match": True,
                    "original_member_filename_in_public_metadata": TARGET_FILENAME in {m["filename"] for m in members},
                    "public_file_metadata": members,
                }
        except (HTTPError, URLError, TimeoutError, ValueError, json.JSONDecodeError) as e:
            attempts.append({"query": query, "result": "ACCESS_OR_SCHEMA_ERROR",
                             "error_type": type(e).__name__, "http_status": getattr(e,"code",None)})
    records = [matches[k] for k in sorted(matches)]
    return {
        "schema": "ttf_Bontrager2021_Zenodo_original_Dryad_public_preservation_index_v01",
        "original_dryad_doi": TARGET_DOI,
        "expected_original_filename": TARGET_FILENAME,
        "status": ("EXACT_DOI_PUBLIC_ZENODO_PRESERVATION_RECORD_LOCATED"
                   if records else "NO_EXACT_DOI_ZENODO_RECORD_DISCOVERED_NOT_ABSENCE_PROOF"),
        "matched_records": records,
        "search_attempts": attempts,
        "absence_in_index_NOT_proof_of_archive_absence": True,
        "fitness_or_experiment_data_read": False,
        "source_study_experiment_overlap_measured": False,
        "Bontrager_44_of_50_frozen_stop_unchanged": True,
    }

if __name__=="__main__":
    print("BONTRAGER_ZENODO_ARCHIVE_PREFLIGHT="+json.dumps(probe(),sort_keys=True))
