"""Pinned original Dryad v3 Timema source-only pair-geometry gate.

No host preferences, mating outcomes, or genetic Fst values are read into models.
The source public file SHA256 is verified before any metadata parsing.
"""
from __future__ import annotations
import hashlib
import json
import urllib.request
import urllib.error

from ttf_crest_timema_source_only_20261115 import audit_bytes

VERSION_ID = 275058
VERSION_NUMBER = 3
PUBLISHED = "2024-01-23"
FILE_CONTRACT = {
    "HostRIDat.csv": {
        "id": 2860173,
        "sha256": "c619bbebbb2a8b84a90da3eafbad2110c0c9fb8a685027bee5c5563ddadee521",
        "size": 5495,
    },
    "CountData.csv": {
        "id": 2860162,
        "sha256": "80e02a00fdf0cdd915c167269c5c8692d92e9fa5e991ef1926c20072994fdef4",
        "size": 1556,
    },
}


def fetch_pinned_public_file(name: str) -> bytes:
    cfg = FILE_CONTRACT[name]
    url = f"https://datadryad.org/api/v2/files/{cfg['id']}/download"
    request = urllib.request.Request(url, headers={"User-Agent": "TTF-source-only/1.0"})
    with urllib.request.urlopen(request, timeout=50) as response:
        raw = response.read(cfg["size"] + 1)
    if len(raw) != cfg["size"] or hashlib.sha256(raw).hexdigest() != cfg["sha256"]:
        raise ValueError(f"published Dryad v3 source byte/digest mismatch: {name}")
    return raw


def main():
    try:
        host = fetch_pinned_public_file("HostRIDat.csv")
        mating = fetch_pinned_public_file("CountData.csv")
    except urllib.error.HTTPError as e:
        if e.code != 401:
            raise
        print(json.dumps({
            "status": "RAW_SOURCE_NOT_ACCESSIBLE_AUTH_REQUIRED_NOT_EVALUABLE",
            "http_status": 401,
            "source": {"doi":"10.5061/dryad.q573n5tpk",
                       "version_id": VERSION_ID, "version_number": VERSION_NUMBER,
                       "published": PUBLISHED, "files": FILE_CONTRACT},
            "pair_counts": None,
            "species_disjoint_holdout": None,
            "no_biological_outcomes_opened": True,
            "not_a_genetic_or_reproductive_isolation_null": True,
        }, sort_keys=True, indent=2))
        return
    result = audit_bytes(host, mating,
                         FILE_CONTRACT["HostRIDat.csv"]["sha256"],
                         FILE_CONTRACT["CountData.csv"]["sha256"])
    result["source"] = {
        "doi": "10.5061/dryad.q573n5tpk",
        "version_id": VERSION_ID,
        "version_number": VERSION_NUMBER,
        "publication_date": PUBLISHED,
        "files": FILE_CONTRACT,
    }
    result["no_raw_source_files_exported"] = True
    print(json.dumps(result, ensure_ascii=False, sort_keys=True, indent=2))


if __name__ == "__main__":
    main()
