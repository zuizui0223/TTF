#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
from pathlib import Path
import urllib.request

RECORDS = (7626793, 7626736)


def fetch(record_id: int) -> dict:
    url = f"https://zenodo.org/api/records/{record_id}"
    request = urllib.request.Request(url, headers={"User-Agent": "TTF Amazon66 mirror probe/0.1", "Accept": "application/json"})
    with urllib.request.urlopen(request, timeout=90) as response:
        payload = json.loads(response.read().decode("utf-8"))
    files = []
    for item in payload.get("files", []):
        files.append({
            "key": item.get("key"),
            "size": item.get("size"),
            "checksum": item.get("checksum"),
            "links": item.get("links", {}),
        })
    return {
        "record_id": record_id,
        "doi": payload.get("doi"),
        "conceptdoi": payload.get("conceptdoi"),
        "title": payload.get("metadata", {}).get("title"),
        "resource_type": payload.get("metadata", {}).get("resource_type"),
        "publication_date": payload.get("metadata", {}).get("publication_date"),
        "files": files,
        "links": payload.get("links", {}),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description="Probe the two Zenodo records associated with the Amazon66 dataset/code lineage without downloading files.")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    rows = [fetch(record_id) for record_id in RECORDS]
    target_names = {"River_islands_sampling_masked.csv", "README.md", "Variable_key.csv", "vcf_files.zip"}
    summary = []
    for row in rows:
        names = {str(item.get("key")) for item in row["files"]}
        summary.append({
            "record_id": row["record_id"],
            "doi": row["doi"],
            "title": row["title"],
            "n_files": len(row["files"]),
            "target_files_present": sorted(target_names & names),
        })
    payload = {
        "schema": "ttf_amazon66_zenodo_probe_v0.1",
        "selection": {
            "file_contents_downloaded": False,
            "empirical_genetic_outcomes_used": False,
            "known_boundary_labels_used": False,
        },
        "records": rows,
        "summary": summary,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(summary, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
