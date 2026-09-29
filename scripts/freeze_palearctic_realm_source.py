#!/usr/bin/env python3
from __future__ import annotations

import argparse
import hashlib
import json
import tempfile
import zipfile
from pathlib import Path


SOURCE_URL = "https://storage.googleapis.com/teow2016/Ecoregions2017.zip"
SOURCE_DOI = "10.1093/biosci/bix014"


def sha256_path(path: Path) -> str:
    h=hashlib.sha256()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda:handle.read(1<<20),b""):
            h.update(chunk)
    return h.hexdigest()


def main() -> int:
    ap=argparse.ArgumentParser()
    ap.add_argument("--zip",type=Path,required=True)
    ap.add_argument("--output",type=Path,required=True)
    args=ap.parse_args()

    archive=args.zip.resolve()
    if not archive.is_file():
        raise FileNotFoundError(archive)

    try:
        import shapefile
    except ImportError as exc:
        raise RuntimeError("pyshp is required to inspect Ecoregions2017") from exc

    with tempfile.TemporaryDirectory(prefix="ttf_ecoregions2017_") as td:
        root=Path(td)
        with zipfile.ZipFile(archive) as zf:
            zf.extractall(root)
        shp=list(root.rglob("Ecoregions2017.shp"))
        if len(shp)!=1:
            raise RuntimeError(f"expected exactly one Ecoregions2017.shp, found {len(shp)}")
        shp=shp[0]
        stem=shp.with_suffix("")
        required=[stem.with_suffix(ext) for ext in (".shp",".shx",".dbf",".prj")]
        missing=[str(p.name) for p in required if not p.is_file()]
        if missing:
            raise RuntimeError(f"Ecoregions2017 vector components missing: {missing}")

        reader=shapefile.Reader(str(shp))
        fields=[row[0] for row in reader.fields[1:]]
        if "REALM" not in fields:
            raise RuntimeError("Ecoregions2017 lacks required REALM field")
        realm_idx=fields.index("REALM")
        realms=[]
        palearctic=0
        for rec in reader.iterRecords():
            value=str(rec[realm_idx])
            realms.append(value)
            if value=="Palearctic":
                palearctic += 1
        if palearctic<=0:
            raise RuntimeError("Ecoregions2017 contains no Palearctic records")

        components={
            p.suffix.lower():{
                "filename":p.name,
                "size_bytes":p.stat().st_size,
                "sha256":sha256_path(p),
            }
            for p in required
        }
        payload={
            "schema":"ttf_palearctic_realm_source_receipt_v0.1",
            "status":"FROZEN_BEFORE_PANEL_MEMBERSHIP",
            "source_url":SOURCE_URL,
            "citation_doi":SOURCE_DOI,
            "archive":{
                "filename":archive.name,
                "size_bytes":archive.stat().st_size,
                "sha256":sha256_path(archive),
            },
            "vector":{
                "dataset":"Ecoregions2017",
                "feature_count":len(reader),
                "realm_field":"REALM",
                "unique_realms":sorted(set(realms)),
                "palearctic_feature_count":palearctic,
                "components":components,
                "ecoregions_vector_sha256":components[".shp"]["sha256"],
            },
            "panel_membership_opened":False,
            "genetic_pair_response_opened":False,
            "interpretation":"This freezes only the external realm geometry source before any authoritative Palearctic species membership is computed."
        }

    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(payload,indent=2,sort_keys=True)+"\n")
    print(json.dumps({
        "status":payload["status"],
        "archive_sha256":payload["archive"]["sha256"],
        "feature_count":payload["vector"]["feature_count"],
        "palearctic_feature_count":payload["vector"]["palearctic_feature_count"],
        "unique_realms":payload["vector"]["unique_realms"],
        "ecoregions_vector_sha256":payload["vector"]["ecoregions_vector_sha256"],
        "panel_membership_opened":False,
    },sort_keys=True))
    return 0


if __name__=="__main__":
    raise SystemExit(main())
