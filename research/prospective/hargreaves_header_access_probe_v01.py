"""Header-only HTTP access preflight for public Hargreaves et al 2020 Dryad data.

No biological response is opened; the first CSV header line is the ONLY content
parsed. This preflight does not invent source-record overlap from dataset pages.
The requested URLs come from Dryad's public Hargreaves dataset file listing.
"""
from __future__ import annotations

import csv
import hashlib
import io
import json
from urllib.error import HTTPError, URLError
from urllib.request import Request, urlopen

FILE_ID = 244459
DATASET_DOI = "10.5061/dryad.0vt4b8gv1"
EXPECTED_FILENAME = "Hargreaves_etal_2020_AmNat_transplant_data.csv"
URLS = (
    f"https://datadryad.org/api/v2/files/{FILE_ID}/download",
    f"https://datadryad.org/downloads/file_stream/{FILE_ID}",
)
MAX_HEADER_BYTES = 32768


def parse_header_bytes(initial: bytes) -> list[str]:
    """Strictly consume only the first header line, never a numerical data row."""
    if len(initial) > MAX_HEADER_BYTES:
        raise ValueError("Header byte cap exceeded")
    if initial.startswith(b"<html") or initial.startswith(b"<!DOCTYPE"):
        raise ValueError("HTML challenge cannot be parsed as CSV")
    line = initial.split(b"\n", 1)[0].rstrip(b"\r")
    if len(line) > 10000:
        raise ValueError("Implausible CSV header length")
    cells = next(csv.reader(io.StringIO(line.decode("utf-8-sig"))))
    if not 2 <= len(cells) <= 500 or len(set(cells)) != len(cells):
        raise ValueError("Malformed CSV column schema")
    return cells


def probe() -> dict:
    attempts = []
    for url in URLS:
        req = Request(
            url,
            headers={
                "User-Agent": "TTF-source-only-preflight/2026.11 (+source-schema)",
                "Accept": "text/csv,application/octet-stream,*/*",
                "Range": f"bytes=0-{MAX_HEADER_BYTES - 1}",
            },
        )
        try:
            with urlopen(req, timeout=35) as response:
                data = response.read(MAX_HEADER_BYTES)
                header = parse_header_bytes(data)
                return {
                    "schema": "ttf_Hargreaves2020_source_only_header_preflight_v01",
                    "dataset_doi": DATASET_DOI,
                    "file_id": FILE_ID,
                    "filename": EXPECTED_FILENAME,
                    "status": "PUBLIC_CSV_HEADER_OBTAINED_METADATA_FIELDS_ONLY",
                    "url_kind": "Dryad API" if "/api/" in url else "Dryad direct",
                    "http_status": int(response.status),
                    "column_count": len(header),
                    "columns": header,
                    "header_sha256": hashlib.sha256(
                        ",".join(header).encode("utf-8")
                    ).hexdigest(),
                    "fitness_rows_read": 0,
                    "source_study_overlap_checked": False,
                    "records_claimed_independent": 0,
                    "attempts": attempts,
                }
        except (HTTPError, URLError, TimeoutError, ValueError, UnicodeError) as exc:
            if isinstance(exc, HTTPError):
                reason = f"HTTP_{exc.code}"
            elif isinstance(exc, (URLError, TimeoutError)):
                reason = "NETWORK_OR_TIMEOUT"
            else:
                reason = "NON_CSV_RESPONSE"
            attempts.append({
                "url_kind": "Dryad API" if "/api/" in url else "Dryad direct",
                "reason": reason,
            })
    return {
        "schema": "ttf_Hargreaves2020_source_only_header_preflight_v01",
        "dataset_doi": DATASET_DOI,
        "file_id": FILE_ID,
        "filename": EXPECTED_FILENAME,
        "status": "SOURCE_ACCESS_GATED_NO_RECORD_LEVEL_OVERLAP_CLAIM",
        "attempts": attempts,
        "fitness_rows_read": 0,
        "source_study_overlap_checked": False,
        "records_claimed_independent": 0,
    }


if __name__ == "__main__":
    print("TTF_HARGREAVES_SOURCE_PREFLIGHT=" + json.dumps(probe(), sort_keys=True))
