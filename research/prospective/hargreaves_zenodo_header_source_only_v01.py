"""Zenodo-preserved Hargreaves source CSV: response-blind FIRST-LINE-ONLY probe.

The exact Dryad DOI must be present as a structured identifier in original
Zenodo metadata, not inferred from a title. Read just the CSV header line
through an official zenodo.org file-content URL. Do not consume any data row.
A successful header preflight does NOT certify original primary-study overlap.
"""
from __future__ import annotations

import csv
import hashlib
import io
import json
from urllib.parse import urlparse
from urllib.request import Request, urlopen

from hargreaves_zenodo_preservation_index_v01 import (
    TARGET_DOI, TARGET_FILENAME, is_original_dryad_doimatch,
)

RECORD_ID = 4409291
RECORD_URL = f"https://zenodo.org/api/records/{RECORD_ID}"
MAX_RECORD_METADATA_BYTES = 250_000
MAX_HEADER_BYTES = 32_768


def original_file_info(record: dict) -> dict:
    if str(record.get("id")) != str(RECORD_ID) or not is_original_dryad_doimatch(record):
        raise ValueError("Record ID and original Dryad DOI must both be exact")
    group = record.get("files") or {}
    entries = group.get("entries") if isinstance(group, dict) else None
    if not isinstance(entries, dict):
        raise ValueError("Expected public Zenodo files.entries metadata")
    matches = [
        v for key, v in entries.items()
        if isinstance(v, dict) and (v.get("key") or key) == TARGET_FILENAME
    ]
    if len(matches) != 1:
        raise ValueError("Original Hargreaves CSV absent or ambiguous")
    f = matches[0]
    link = (f.get("links") or {}).get("content")
    if not isinstance(link, str):
        raise ValueError("No public metadata-provided official content URL")
    u = urlparse(link)
    if (u.scheme != "https" or u.hostname != "zenodo.org"
        or u.username is not None or u.password is not None
        or not u.path.startswith(f"/api/records/{RECORD_ID}/files/")
        or not u.path.endswith("/content")):
        raise ValueError("Source-file content URL is not verified official Zenodo API")
    size = f.get("size")
    if not isinstance(size, int) or not 1_000_000 <= size <= 6_000_000:
        raise ValueError("Original CSV file size outside expected published metadata range")
    return {"content_url": link, "size": size,
            "checksum_from_published_metadata": f.get("checksum")}


def parse_header(line: bytes) -> list[str]:
    if not line or len(line) > MAX_HEADER_BYTES or b"\n" not in line:
        raise ValueError("No complete bounded CSV header")
    header = next(csv.reader(io.StringIO(line.decode("utf-8-sig")),
                             delimiter=",", strict=True))
    if not 5 <= len(header) <= 500 or not all(x.strip() for x in header):
        raise ValueError("Implausible original source CSV header")
    if len(set(header)) != len(header):
        raise ValueError("Duplicate source CSV columns")
    return header


def probe() -> dict:
    result = {
        "schema": "ttf_Hargreaves2020_Zenodo_exact_original_source_header_only_v01",
        "original_dryad_doi": TARGET_DOI,
        "zenodo_preservation_record_id": RECORD_ID,
        "original_csv_filename": TARGET_FILENAME,
        "source_study_overlap_measured": False,
        "independent_new_studies_certified": 0,
        "fitness_rows_read": 0,
        "response_values_modeled": False,
    }
    try:
        req = Request(RECORD_URL, headers={"Accept": "application/json",
            "User-Agent": "TTF-Hargreaves-official-preservation-header-only/2026-11-15"})
        with urlopen(req, timeout=30) as response:
            metadata_bytes = response.read(MAX_RECORD_METADATA_BYTES + 1)
        if len(metadata_bytes) > MAX_RECORD_METADATA_BYTES:
            raise ValueError("Record metadata exceeds cap")
        record = json.loads(metadata_bytes)
        file = original_file_info(record)
        req = Request(file["content_url"], headers={
            "Accept": "text/csv,application/octet-stream",
            "User-Agent": "TTF-source-header-firstline-only/2026-11-15"})
        with urlopen(req, timeout=30) as response:
            header_line = response.readline(MAX_HEADER_BYTES + 1)
        header = parse_header(header_line)
        result.update({
            "status": "OFFICIAL_PRESERVED_ORIGINAL_CSV_HEADER_OBTAINED_SOURCE_ONLY",
            "source_csv_columns": header,
            "source_csv_column_count": len(header),
            "source_csv_header_SHA256": hashlib.sha256(header_line).hexdigest(),
            "zenodo_file_size_bytes": file["size"],
            "zenodo_file_checksum_metadata": file["checksum_from_published_metadata"],
            "original_csv_data_rows_read": 0,
            "fitness_response_unopened": True,
        })
    except Exception as exc:
        result.update({
            "status": "OFFICIAL_ZENODO_PRESERVED_CSV_HEADER_ACCESS_NOT_EVALUABLE",
            "error_type": type(exc).__name__,
            "error_http_code": getattr(exc, "code", None),
            "error_message": str(exc)[:220] if not hasattr(exc, "code") else None,
        })
    return result


if __name__ == "__main__":
    print("HARGREAVES_ZENODO_HEADER=" + json.dumps(probe(), sort_keys=True))
