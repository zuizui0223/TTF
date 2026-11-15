"""Inspect PUBLIC Dryad version/file manifests, metadata ONLY; no data files fetched.

Source dataset: doi:10.5061/dryad.q573n5tpk.
This preflight must not read host-choice, mating counts, or genetic outcomes.
"""
from __future__ import annotations
import json
import urllib.parse
import urllib.request

DOI = "doi:10.5061/dryad.q573n5tpk"
BASE = "https://datadryad.org/api/v2"
ENCODED = urllib.parse.quote(DOI, safe="")


def fetch_json(url):
    request = urllib.request.Request(
        url,
        headers={"Accept": "application/json", "User-Agent": "TTF-research-source-preflight/1.0"},
    )
    with urllib.request.urlopen(request, timeout=40) as response:
        if response.status != 200:
            raise ValueError(f"API status {response.status}")
        raw = response.read(3_000_000)
    return json.loads(raw)


def nested_list(document, keyword):
    embedded = document.get("_embedded", {})
    for key, value in embedded.items():
        if keyword.lower() in key.lower() and isinstance(value, list):
            return value, key
    return [], ""


def project_link(link):
    if isinstance(link, str): return link
    if isinstance(link, dict): return link.get("href", "")
    return ""


def safe_version(v):
    links = v.get("_links", {})
    href = project_link(links.get("self", ""))
    id_from_self = int(href.rsplit("/",1)[-1]) if href.rsplit("/",1)[-1].isdigit() else None
    return {
        "id": v.get("id") or id_from_self,
        "versionNumber": v.get("versionNumber"),
        "publicationDate": v.get("publicationDate"),
        "lastModificationDate": v.get("lastModificationDate"),
        "date": v.get("date"),
        "self_href": project_link(links.get("self", "")),
    }


def safe_file(f):
    links = f.get("_links", {})
    return {
        "id": f.get("id") or (int(project_link(links.get("self", "")).rsplit("/",1)[-1]) if project_link(links.get("self", "")).rsplit("/",1)[-1].isdigit() else None),
        "fileName": f.get("fileName") or f.get("path"),
        "path": f.get("path"),
        "size": f.get("size"),
        "digest": f.get("digest"),
        "digestType": f.get("digestType"),
        "md5": f.get("md5"),
        "lastModificationDate": f.get("lastModificationDate"),
        "self_href": project_link(links.get("self", "")),
        "download_href": project_link(links.get("stash:download", links.get("download", ""))),
    }


def main():
    versions_url = f"{BASE}/datasets/{ENCODED}/versions"
    a = fetch_json(versions_url)
    versions, version_collection = nested_list(a, "versions")
    if not versions:
        raise ValueError(f"versions not exposed in expected HAL embedded collection: keys={list(a.get('_embedded',{}))}")
    version_rows = [safe_version(x) for x in versions]
    candidates = [x for x in version_rows if isinstance(x["id"], int)]
    if not candidates:
        raise ValueError(f"no numeric version ID: {version_rows}")
    # Deterministic latest-number metadata for the *discovery* run.
    versions_2024 = [x for x in candidates if x["versionNumber"] == 3 and x["id"] == 275058
                     and x["publicationDate"] == "2024-01-23"]
    if len(versions_2024) != 1:
        raise ValueError("published 2024 Dryad version 3 ID 275058 not uniquely available")
    v = versions_2024[0]
    files_url = f"{BASE}/versions/{v['id']}/files"
    b = fetch_json(files_url)
    files, files_collection = nested_list(b, "files")
    target_names = {"HostRIDat.csv", "CountData.csv", "README.md"}
    result = {
        "status": "DRYAD_METADATA_ONLY_NO_DATA_FILE_ACCESSED",
        "doi": DOI, "versions_url": versions_url,
        "version_collection": version_collection,
        "versions": version_rows,
        "selected_version_metadata_only": v,
        "files_url": files_url,
        "file_collection": files_collection,
        "files_total_page": len(files),
        "target_file_metadata": [safe_file(x) for x in files if (x.get("fileName") or x.get("path")) in target_names],
        "available_filenames": [x.get("fileName") or x.get("path") for x in files],
        "no_file_contents_fetched": True,
        "no_outcomes_opened": True,
    }
    if {x["fileName"] for x in result["target_file_metadata"]} != target_names:
        summary = {"file_embedded_keys": list(b.get("_embedded",{})),
                   "first_file_keys": list(files[0]) if files else [],
                   "first_file_preview": {k: str(v)[:180] for k,v in (files[0] if files else {}).items() if k not in ("description",)},
                   "file_count": len(files), "root_keys": list(b)}
        raise ValueError("Dryad target names unresolved (metadata schema diagnosis): " + json.dumps(summary))
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
