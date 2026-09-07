#!/usr/bin/env python3
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import shutil
import urllib.parse
import urllib.request

DOI = "10.5061/dryad.rxwdbrvc1"
DOI_ID = f"doi:{DOI}"
API = "https://datadryad.org/api/v2"
LANDING = "https://datadryad.org/dataset/doi:10.5061/dryad.rxwdbrvc1"
REQUIRED_SMALL = (
    "README.md",
    "River_islands_sampling_masked.csv",
    "Variable_key.csv",
    "trait.database.formatted.final.csv",
)
GENETIC_ARCHIVE = "vcf_files.zip"


def request_json(url: str) -> dict:
    req = urllib.request.Request(
        url,
        headers={"User-Agent": "TTF Amazon66 metadata preflight/0.1", "Accept": "application/json"},
    )
    with urllib.request.urlopen(req, timeout=90) as response:
        return json.loads(response.read().decode("utf-8"))


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def collect_dicts(value):
    if isinstance(value, dict):
        yield value
        for child in value.values():
            yield from collect_dicts(child)
    elif isinstance(value, list):
        for child in value:
            yield from collect_dicts(child)


def hrefs_from_links(value: dict) -> list[tuple[str, str]]:
    out: list[tuple[str, str]] = []
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


def select_latest_version(dataset: dict) -> tuple[int, dict]:
    # Prefer a singular version relation if Dryad exposes one.
    for key, href in hrefs_from_links(dataset):
        if "version" in key.lower() and not href.rstrip("/").endswith("versions"):
            candidate = request_json(absolute(href))
            if isinstance(candidate.get("id"), int):
                return int(candidate["id"]), candidate

    encoded = urllib.parse.quote(DOI_ID, safe="")
    versions_url = f"{API}/datasets/{encoded}/versions"
    payload = request_json(versions_url)
    candidates = []
    for item in collect_dicts(payload):
        if isinstance(item.get("id"), int) and (
            "versionNumber" in item or "versionStatus" in item or "version" in item
        ):
            candidates.append(item)
    if not candidates:
        raise RuntimeError(f"could not resolve a Dryad version from {versions_url}")
    candidates.sort(key=lambda item: (int(item.get("versionNumber", 0)), int(item["id"])))
    chosen = candidates[-1]
    return int(chosen["id"]), chosen


def file_inventory(files_payload: dict) -> list[dict]:
    seen = set()
    rows: list[dict] = []
    for item in collect_dicts(files_payload):
        path = item.get("path") or item.get("name") or item.get("filename")
        if not isinstance(path, str):
            continue
        basename = Path(path).name
        if not basename:
            continue
        size = item.get("size")
        if not isinstance(size, int):
            continue
        ident = item.get("id")
        key = (ident, path, size)
        if key in seen:
            continue
        seen.add(key)
        download_candidates = []
        for field in ("downloadURL", "downloadUrl", "downloadLink", "download_url"):
            if isinstance(item.get(field), str):
                download_candidates.append(item[field])
        for rel, href in hrefs_from_links(item):
            if "download" in rel.lower() or "content" in rel.lower():
                download_candidates.append(href)
        rows.append(
            {
                "id": ident,
                "path": path,
                "name": basename,
                "size": size,
                "digest": item.get("digest"),
                "digest_type": item.get("digestType") or item.get("digest_type"),
                "mime_type": item.get("mimeType") or item.get("mime_type"),
                "download_candidates": list(dict.fromkeys(map(str, download_candidates))),
            }
        )
    rows.sort(key=lambda row: row["name"])
    return rows


def download_candidates(row: dict) -> list[str]:
    candidates = [absolute(url) for url in row.get("download_candidates", [])]
    if row.get("id") is not None:
        candidates.append(f"https://datadryad.org/stash/downloads/file_stream/{row['id']}")
    return list(dict.fromkeys(candidates))


def try_download(row: dict, destination: Path) -> dict:
    attempts = []
    for url in download_candidates(row):
        try:
            req = urllib.request.Request(
                url,
                headers={"User-Agent": "TTF Amazon66 metadata preflight/0.1"},
            )
            with urllib.request.urlopen(req, timeout=120) as response, destination.open("wb") as out:
                shutil.copyfileobj(response, out)
            if destination.stat().st_size != int(row["size"]):
                raise RuntimeError(
                    f"size mismatch: expected {row['size']}, got {destination.stat().st_size}"
                )
            return {
                "downloaded": True,
                "url": url,
                "bytes": destination.stat().st_size,
                "sha256": sha256(destination),
                "attempts": attempts,
            }
        except Exception as exc:  # retain exact failure and try only equivalent Dryad links
            attempts.append({"url": url, "error": f"{type(exc).__name__}: {exc}"})
            if destination.exists():
                destination.unlink()
    return {"downloaded": False, "url": None, "bytes": None, "sha256": None, "attempts": attempts}


def main() -> int:
    parser = argparse.ArgumentParser(description="Freeze metadata-only preflight for the Johnson et al. Amazonian 66-species Dryad dataset.")
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    output = args.output_dir
    output.mkdir(parents=True, exist_ok=True)

    encoded = urllib.parse.quote(DOI_ID, safe="")
    dataset_url = f"{API}/datasets/{encoded}"
    dataset = request_json(dataset_url)
    version_id, version = select_latest_version(dataset)
    files_url = f"{API}/versions/{version_id}/files?per_page=100"
    files_payload = request_json(files_url)
    inventory = file_inventory(files_payload)
    by_name: dict[str, list[dict]] = {}
    for row in inventory:
        by_name.setdefault(row["name"], []).append(row)

    required_status = {}
    for name in REQUIRED_SMALL:
        matches = by_name.get(name, [])
        if len(matches) != 1:
            required_status[name] = {
                "present_exactly_once": False,
                "matches": matches,
                "download": None,
            }
            continue
        row = matches[0]
        destination = output / name
        required_status[name] = {
            "present_exactly_once": True,
            "file": row,
            "download": try_download(row, destination),
        }

    vcf_matches = by_name.get(GENETIC_ARCHIVE, [])
    # Deliberately do not download genetic outcomes during geometry preflight.
    vcf_status = {
        "present_exactly_once": len(vcf_matches) == 1,
        "matches": vcf_matches,
        "download_attempted": False,
    }

    ready = all(
        status.get("present_exactly_once")
        and status.get("download", {}).get("downloaded")
        for status in required_status.values()
    ) and vcf_status["present_exactly_once"]

    manifest = {
        "schema": "ttf_amazon66_metadata_preflight_v0.1",
        "source": {
            "dataset_doi": DOI,
            "landing_page": LANDING,
            "dataset_api": dataset_url,
            "version_id": version_id,
            "files_api": files_url,
        },
        "selection": {
            "empirical_genetic_outcomes_used": False,
            "known_boundary_labels_used": False,
            "vcf_download_attempted": False,
            "species_filtering_on_genetic_results": False,
        },
        "dataset_metadata": {
            "title": dataset.get("title"),
            "version_number": dataset.get("versionNumber") or version.get("versionNumber"),
            "storage_size": dataset.get("storageSize"),
            "license": dataset.get("license"),
        },
        "file_inventory": inventory,
        "required_small_files": required_status,
        "genetic_archive": vcf_status,
        "ready_for_geometry_audit": bool(ready),
    }
    (output / "preflight_manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps({
        "version_id": version_id,
        "n_files": len(inventory),
        "required": {
            name: {
                "present": status.get("present_exactly_once"),
                "downloaded": None if status.get("download") is None else status["download"].get("downloaded"),
            }
            for name, status in required_status.items()
        },
        "vcf_present": vcf_status["present_exactly_once"],
        "vcf_download_attempted": False,
        "ready_for_geometry_audit": ready,
    }, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
