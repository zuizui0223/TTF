#!/usr/bin/env python3
from __future__ import annotations

import argparse
import hashlib
import http.cookiejar
import json
from pathlib import Path
import urllib.error
import urllib.request

LANDING = "https://datadryad.org/dataset/doi:10.5061/dryad.m7rc3"
FILE_ID = 46189
EXPECTED_SIZE = 31162
EXPECTED_MD5 = "e2ca197cda27d6f4c10b02f5a20f41da"
CANDIDATES = (
    f"https://datadryad.org/stash/downloads/file_stream/{FILE_ID}",
    f"https://datadryad.org/stash/downloads/file_stream/{FILE_ID}?download=1",
)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    cookies = http.cookiejar.CookieJar()
    opener = urllib.request.build_opener(urllib.request.HTTPCookieProcessor(cookies))
    base_headers = {
        "User-Agent": "Mozilla/5.0 (X11; Linux x86_64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/140.0 Safari/537.36",
        "Accept-Language": "en-US,en;q=0.9",
    }
    landing_status = None
    landing_error = None
    try:
        req = urllib.request.Request(LANDING, headers={**base_headers, "Accept": "text/html,application/xhtml+xml"})
        with opener.open(req, timeout=90) as response:
            landing_status = int(response.status)
            response.read(4096)
    except Exception as exc:
        landing_error = f"{type(exc).__name__}: {exc}"

    attempts = []
    success = None
    for url in CANDIDATES:
        try:
            req = urllib.request.Request(
                url,
                headers={
                    **base_headers,
                    "Accept": "application/zip,application/octet-stream,*/*",
                    "Referer": LANDING,
                    "Sec-Fetch-Dest": "document",
                    "Sec-Fetch-Mode": "navigate",
                    "Sec-Fetch-Site": "same-origin",
                },
            )
            with opener.open(req, timeout=120) as response:
                body = response.read()
                status = int(response.status)
                content_type = response.headers.get("Content-Type")
            md5 = hashlib.md5(body).hexdigest()
            attempts.append({"url": url, "status": status, "content_type": content_type, "bytes": len(body), "md5": md5})
            if len(body) == EXPECTED_SIZE and md5 == EXPECTED_MD5:
                success = attempts[-1]
                break
        except Exception as exc:
            code = exc.code if isinstance(exc, urllib.error.HTTPError) else None
            attempts.append({"url": url, "status": code, "error": f"{type(exc).__name__}: {exc}"})

    payload = {
        "schema": "ttf_queensland33_session_download_probe_v0.1",
        "landing": {"url": LANDING, "status": landing_status, "error": landing_error, "cookie_names": sorted(cookie.name for cookie in cookies)},
        "target": {"file_id": FILE_ID, "expected_size": EXPECTED_SIZE, "expected_md5": EXPECTED_MD5},
        "attempts": attempts,
        "matched_frozen_file": success is not None,
        "content_persisted": False,
        "selection": {
            "sequence_contents_reported": False,
            "genetic_distances_computed": False,
            "known_boundary_labels_used": False,
        },
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    print(json.dumps(payload, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
