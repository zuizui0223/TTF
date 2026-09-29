#!/usr/bin/env python3
from __future__ import annotations

import argparse
import hashlib
import json
import tempfile
import zipfile
from pathlib import Path


def sha256_path(path: Path) -> str:
    h=hashlib.sha256()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda:handle.read(1<<20),b""):
            h.update(chunk)
    return h.hexdigest()


def main() -> int:
    ap=argparse.ArgumentParser()
    ap.add_argument("--zip",type=Path,required=True)
    ap.add_argument("--source-receipt",type=Path,required=True)
    ap.add_argument("--output-geojson",type=Path,required=True)
    ap.add_argument("--output-receipt",type=Path,required=True)
    args=ap.parse_args()

    source=json.loads(args.source_receipt.read_text())
    if source.get("schema")!="ttf_palearctic_realm_source_receipt_v0.1":
        raise RuntimeError("unexpected Palearctic realm source receipt")
    if source.get("status")!="FROZEN_BEFORE_PANEL_MEMBERSHIP":
        raise RuntimeError("realm source receipt is not pre-membership frozen")
    if sha256_path(args.zip)!=source["archive"]["sha256"]:
        raise RuntimeError("Ecoregions2017 archive SHA256 differs from frozen receipt")

    try:
        import shapefile
    except ImportError as exc:
        raise RuntimeError("pyshp is required") from exc

    with tempfile.TemporaryDirectory(prefix="ttf_palearctic_subset_") as td:
        root=Path(td)
        with zipfile.ZipFile(args.zip) as zf:
            zf.extractall(root)
        shp=list(root.rglob("Ecoregions2017.shp"))
        if len(shp)!=1:
            raise RuntimeError(f"expected one Ecoregions2017.shp, found {len(shp)}")
        reader=shapefile.Reader(str(shp[0]),encoding="latin1")
        fields=[row[0] for row in reader.fields[1:]]
        for required in ("REALM","ECO_NAME"):
            if required not in fields:
                raise RuntimeError(f"required Ecoregions2017 field missing: {required}")
        realm_idx=fields.index("REALM")
        eco_name_idx=fields.index("ECO_NAME")
        eco_id_idx=fields.index("ECO_ID") if "ECO_ID" in fields else None

        features=[]
        for record,shape in zip(reader.iterRecords(),reader.iterShapes()):
            if str(record[realm_idx])!="Palearctic":
                continue
            properties={
                "REALM":"Palearctic",
                "ECO_NAME":str(record[eco_name_idx]),
            }
            if eco_id_idx is not None:
                properties["ECO_ID"]=record[eco_id_idx]
            features.append({
                "type":"Feature",
                "properties":properties,
                "geometry":shape.__geo_interface__,
            })

    if len(features)!=int(source["vector"]["palearctic_feature_count"]):
        raise RuntimeError(
            "Palearctic feature count differs from frozen source receipt"
        )
    features.sort(key=lambda f:(str(f["properties"].get("ECO_ID","")),f["properties"]["ECO_NAME"]))
    payload={
        "type":"FeatureCollection",
        "name":"Ecoregions2017_Palearctic",
        "source_archive_sha256":source["archive"]["sha256"],
        "features":features,
    }
    args.output_geojson.parent.mkdir(parents=True,exist_ok=True)
    args.output_geojson.write_text(
        json.dumps(payload,ensure_ascii=False,sort_keys=True,separators=(",",":"))+"\n"
    )

    receipt={
        "schema":"ttf_palearctic_realm_subset_receipt_v0.1",
        "status":"DERIVED_BEFORE_PANEL_MEMBERSHIP",
        "source_receipt_sha256":sha256_path(args.source_receipt),
        "source_archive_sha256":source["archive"]["sha256"],
        "realm_field":"REALM",
        "realm_value":"Palearctic",
        "feature_count":len(features),
        "geojson_sha256":sha256_path(args.output_geojson),
        "geojson_size_bytes":args.output_geojson.stat().st_size,
        "panel_membership_opened":False,
        "genetic_pair_response_opened":False,
    }
    args.output_receipt.parent.mkdir(parents=True,exist_ok=True)
    args.output_receipt.write_text(json.dumps(receipt,indent=2,sort_keys=True)+"\n")
    print(json.dumps(receipt,sort_keys=True))
    return 0


if __name__=="__main__":
    raise SystemExit(main())
