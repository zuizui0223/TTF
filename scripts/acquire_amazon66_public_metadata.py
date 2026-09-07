#!/usr/bin/env python3
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import shutil
import urllib.request

VERSION_ID = 221232
FILES = {
    "README.md": {
        "file_id": 2113630,
        "bytes": 3356,
        "sha256": "1edc3e139c0dd01491cebf14f7fb40d95279f1b9cf3772ffddef8319fb1850ce",
    },
    "River_islands_sampling_masked.csv": {
        "file_id": 2113627,
        "bytes": 17119,
        "sha256": "7c2bf931946bc63c413d1dff5255d667abbe018a0b15daf161510dcf35363b87",
    },
    "Variable_key.csv": {
        "file_id": 2113628,
        "bytes": 4309,
        "sha256": "7e833d0fe5e8fc0a91ccf4b71b11a0adac4bfa4516cbaf8a1c2044f57a1c83b7",
    },
}
VCF_FROZEN = {
    "name": "vcf_files.zip",
    "file_id": 2113612,
    "bytes": 32478348,
    "sha256": "8f20ef2343a1fcd966745b9bd1fc3353147ca6770d0ed51a7b6faa38ee56e042",
    "download_attempted": False,
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def download_public(file_id: int, destination: Path) -> str:
    url = f"https://datadryad.org/stash/downloads/file_stream/{file_id}"
    request = urllib.request.Request(
        url,
        headers={"User-Agent": "TTF Amazon66 geometry metadata/0.1"},
    )
    with urllib.request.urlopen(request, timeout=120) as response, destination.open("wb") as out:
        shutil.copyfileobj(response, out)
    return url


def main() -> int:
    parser = argparse.ArgumentParser(description="Acquire only frozen non-outcome metadata for Amazon66 geometry qualification.")
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    output = args.output_dir
    output.mkdir(parents=True, exist_ok=True)

    downloaded = []
    errors = []
    for name, expected in FILES.items():
        destination = output / name
        try:
            url = download_public(int(expected["file_id"]), destination)
            observed_bytes = destination.stat().st_size
            observed_sha = sha256(destination)
            if observed_bytes != int(expected["bytes"]):
                raise RuntimeError(f"size mismatch for {name}: {observed_bytes} != {expected['bytes']}")
            if observed_sha != str(expected["sha256"]):
                raise RuntimeError(f"sha256 mismatch for {name}: {observed_sha}")
            downloaded.append({"name": name, "url": url, **expected})
        except Exception as exc:
            if destination.exists():
                destination.unlink()
            errors.append({"name": name, "error": f"{type(exc).__name__}: {exc}"})

    ready = len(downloaded) == len(FILES) and not errors
    manifest = {
        "schema": "ttf_amazon66_public_geometry_metadata_v0.1",
        "source": {
            "dryad_doi": "10.5061/dryad.rxwdbrvc1",
            "dryad_version_id": VERSION_ID,
            "dryad_version_number": 16,
        },
        "selection": {
            "empirical_genetic_outcomes_used": False,
            "published_genetic_summary_tables_opened_for_selection": False,
            "known_boundary_labels_used": False,
            "vcf_download_attempted": False,
        },
        "downloaded": downloaded,
        "errors": errors,
        "genetic_archive_frozen_but_unopened": VCF_FROZEN,
        "ready_for_geometry_audit": ready,
    }
    (output / "geometry_metadata_manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps({"ready_for_geometry_audit": ready, "downloaded": [row["name"] for row in downloaded], "errors": errors}, sort_keys=True))
    if not ready:
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
