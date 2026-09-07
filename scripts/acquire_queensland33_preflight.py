#!/usr/bin/env python3
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import re
import tempfile
import urllib.parse
import urllib.request

DOI = "10.5061/dryad.m7rc3"
DOI_ID = f"doi:{DOI}"
API = "https://datadryad.org/api/v2"
LANDING = "https://datadryad.org/dataset/doi:10.5061/dryad.m7rc3"
TARGETS = (
    "Arlequin input files.zip",
    "MTML-msbayes input files.zip",
    "TCS input files.zip",
)


def request_json(url: str) -> dict:
    req = urllib.request.Request(
        url,
        headers={"User-Agent": "TTF Queensland33 source preflight/0.1", "Accept": "application/json"},
    )
    with urllib.request.urlopen(req, timeout=90) as response:
        return json.loads(response.read().decode("utf-8"))


def collect_dicts(value):
    if isinstance(value, dict):
        yield value
        for child in value.values():
            yield from collect_dicts(child)
    elif isinstance(value, list):
        for child in value:
            yield from collect_dicts(child)


def hrefs(value: dict) -> list[tuple[str, str]]:
    out = []
    links = value.get("_links", {}) if isinstance(value, dict) else {}
    if isinstance(links, dict):
        for key, item in links.items():
            if isinstance(item, dict) and isinstance(item.get("href"), str):
                out.append((str(key), item["href"]))
            elif isinstance(item, list):
                for entry in item:
                    if isinstance(entry, dict) and isinstance(entry.get("href"), str):
                        out.append((str(key), entry["href"]))
    return out


def absolute(url_or_path: str) -> str:
    return urllib.parse.urljoin("https://datadryad.org", url_or_path)


def version_id_from_href(href: str) -> int | None:
    match = re.search(r"/versions/(\d+)(?:$|[/?#])", href)
    return None if match is None else int(match.group(1))


def select_version(dataset: dict) -> tuple[int, dict]:
    for key, link in hrefs(dataset):
        if "version" in key.lower() and not link.rstrip("/").endswith("versions"):
            version_id = version_id_from_href(link)
            payload = request_json(absolute(link))
            if version_id is not None:
                return version_id, payload
            if isinstance(payload.get("id"), int):
                return int(payload["id"]), payload
    encoded = urllib.parse.quote(DOI_ID, safe="")
    versions_url = f"{API}/datasets/{encoded}/versions"
    payload = request_json(versions_url)
    candidates: list[int] = []
    for item in collect_dicts(payload):
        if isinstance(item.get("id"), int):
            candidates.append(int(item["id"]))
        for _key, link in hrefs(item):
            parsed = version_id_from_href(link)
            if parsed is not None:
                candidates.append(parsed)
    if not candidates:
        raise RuntimeError("Dryad version resource could not be resolved")
    version_id = max(candidates)
    return version_id, request_json(f"{API}/versions/{version_id}")


def inventory(payload: dict) -> list[dict]:
    rows = []
    seen = set()
    for item in collect_dicts(payload):
        path = item.get("path") or item.get("name") or item.get("filename")
        size = item.get("size")
        if not isinstance(path, str) or not isinstance(size, int):
            continue
        name = Path(path).name
        key = (name, size, item.get("digest"))
        if key in seen:
            continue
        seen.add(key)
        urls = []
        for field in ("downloadURL", "downloadUrl", "downloadLink", "download_url"):
            if isinstance(item.get(field), str):
                urls.append(absolute(item[field]))
        for rel, link in hrefs(item):
            if "download" in rel.lower() or "content" in rel.lower():
                urls.append(absolute(link))
        # Current Dryad file objects frequently expose the file resource only in
        # the download href, rather than an integer id field.
        file_id = item.get("id") if isinstance(item.get("id"), int) else None
        if file_id is None:
            for url in urls:
                match = re.search(r"/files/(\d+)(?:/|$)", url)
                if match:
                    file_id = int(match.group(1))
                    break
        rows.append({
            "name": name,
            "path": path,
            "size": size,
            "digest": item.get("digest"),
            "digest_type": item.get("digestType") or item.get("digest_type"),
            "file_id": file_id,
            "download_candidates": list(dict.fromkeys(urls)),
        })
    return sorted(rows, key=lambda row: row["name"])


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def probe_download(row: dict) -> dict:
    candidates = list(row.get("download_candidates", []))
    if row.get("file_id") is not None:
        candidates.append(f"https://datadryad.org/stash/downloads/file_stream/{row['file_id']}")
    attempts = []
    with tempfile.TemporaryDirectory(prefix="ttf-qld33-") as tmp:
        destination = Path(tmp) / row["name"]
        for url in dict.fromkeys(candidates):
            try:
                req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0 TTF Queensland33/0.1"})
                with urllib.request.urlopen(req, timeout=120) as response, destination.open("wb") as out:
                    while True:
                        block = response.read(1024 * 1024)
                        if not block:
                            break
                        out.write(block)
                size = destination.stat().st_size
                if size != int(row["size"]):
                    raise RuntimeError(f"size mismatch expected={row['size']} actual={size}")
                return {
                    "downloaded": True,
                    "url": url,
                    "bytes": size,
                    "sha256": sha256(destination),
                    "content_persisted": False,
                    "attempts": attempts,
                }
            except Exception as exc:
                attempts.append({"url": url, "error": f"{type(exc).__name__}: {exc}"})
                if destination.exists():
                    destination.unlink()
    return {
        "downloaded": False,
        "url": None,
        "bytes": None,
        "sha256": None,
        "content_persisted": False,
        "attempts": attempts,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description="Outcome-blind acquisition preflight for Page & Hughes Queensland33 Dryad data.")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    encoded = urllib.parse.quote(DOI_ID, safe="")
    dataset_url = f"{API}/datasets/{encoded}"
    dataset = request_json(dataset_url)
    version_id, version = select_version(dataset)
    files_url = f"{API}/versions/{version_id}/files?per_page=100"
    rows = inventory(request_json(files_url))
    by_name = {row["name"]: row for row in rows}

    target_status = {}
    for name in TARGETS:
        row = by_name.get(name)
        target_status[name] = {
            "present": row is not None,
            "file": row,
            "download_probe": None if row is None else probe_download(row),
        }

    payload = {
        "schema": "ttf_queensland33_source_preflight_v0.1",
        "source": {
            "dataset_doi": DOI,
            "landing_page": LANDING,
            "dataset_api": dataset_url,
            "version_id": version_id,
            "version_number": version.get("versionNumber"),
            "files_api": files_url,
            "dataset_title": dataset.get("title"),
        },
        "selection": {
            "sequence_contents_reported": False,
            "genetic_distances_computed": False,
            "published_divergence_used_for_species_filtering": False,
            "known_Mary_Brisbane_boundary_used_for_filtering": False,
            "downloaded_files_persisted": False,
        },
        "inventory": rows,
        "targets": target_status,
        "ready_for_header_only_parser": all(
            status["present"] and status["download_probe"] and status["download_probe"]["downloaded"]
            for status in target_status.values()
        ),
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps({
        "version_id": version_id,
        "n_files": len(rows),
        "ready_for_header_only_parser": payload["ready_for_header_only_parser"],
        "targets": {
            name: {
                "present": status["present"],
                "downloaded": bool(status["download_probe"] and status["download_probe"]["downloaded"]),
                "size": None if status["file"] is None else status["file"]["size"],
                "file_id": None if status["file"] is None else status["file"]["file_id"],
            }
            for name, status in target_status.items()
        },
    }, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
