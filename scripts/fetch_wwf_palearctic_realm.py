#!/usr/bin/env python3
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import urllib.parse
import urllib.request


SERVICE = (
    "https://data-gis.unep-wcmc.org/server/rest/services/"
    "WWF_004_GeneralisedBiogeographicRealms2004/FeatureServer/0/query"
)


def canonical_query_url() -> str:
    params = {
        "where": "realmcode='PA'",
        "outFields": "realm,realmcode",
        "returnGeometry": "true",
        "outSR": "4326",
        "f": "geojson",
    }
    return SERVICE + "?" + urllib.parse.urlencode(params)


def fetch_bytes(url: str, *, timeout: int = 120) -> bytes:
    request = urllib.request.Request(
        url,
        headers={
            "User-Agent": "TTF-Q-Palearctic-response-blind-realm-fetch/0.1",
            "Accept": "application/geo+json,application/json",
        },
    )
    with urllib.request.urlopen(request, timeout=int(timeout)) as response:
        return response.read()


def validate_geojson(raw: bytes) -> dict:
    payload = json.loads(raw.decode("utf-8"))
    if payload.get("type") != "FeatureCollection":
        raise RuntimeError("UNEP-WCMC realm query did not return a FeatureCollection")
    features = payload.get("features")
    if not isinstance(features, list) or not features:
        raise RuntimeError("UNEP-WCMC realm query returned no features")
    for feature in features:
        props = feature.get("properties") or {}
        code = str(props.get("realmcode", props.get("REALMCODE", ""))).strip()
        name = str(props.get("realm", props.get("REALM", ""))).strip()
        if code != "PA" and name.casefold() != "palearctic":
            raise RuntimeError("realm query returned a non-Palearctic feature")
        geom = feature.get("geometry") or {}
        if geom.get("type") not in {"Polygon", "MultiPolygon"}:
            raise RuntimeError("realm query returned unsupported geometry")
    return payload


def main() -> int:
    ap=argparse.ArgumentParser()
    ap.add_argument("--output-dir",type=Path,required=True)
    args=ap.parse_args()

    url=canonical_query_url()
    raw=fetch_bytes(url)
    payload=validate_geojson(raw)
    out=args.output_dir
    out.mkdir(parents=True,exist_ok=True)
    geojson_path=out/"wwf_palearctic_realm_2004.geojson"
    receipt_path=out/"wwf_palearctic_realm_2004_receipt.json"
    geojson_path.write_bytes(raw)
    receipt={
        "schema":"ttf_wwf_palearctic_realm_fetch_v0.1",
        "status":"FETCHED_AUTHORITATIVE_WWF_PALEARCTIC_GEOJSON",
        "service":SERVICE,
        "query_url":url,
        "realm":"Palearctic",
        "realmcode":"PA",
        "outSR":4326,
        "feature_count":len(payload["features"]),
        "geojson_sha256":hashlib.sha256(raw).hexdigest(),
        "scientific_response_used":False,
    }
    receipt_path.write_text(json.dumps(receipt,indent=2,sort_keys=True)+"\n")
    print(json.dumps({
        "status":receipt["status"],
        "feature_count":receipt["feature_count"],
        "geojson_sha256":receipt["geojson_sha256"],
    },sort_keys=True))
    return 0


if __name__=="__main__":
    raise SystemExit(main())
