#!/usr/bin/env python3
from __future__ import annotations

import argparse
import hashlib
import json
import tempfile
import zipfile
from pathlib import Path

import shapefile


def git_blob_sha1(path: Path) -> str:
    data = Path(path).read_bytes()
    header = f"blob {len(data)}\0".encode("ascii")
    return hashlib.sha1(header + data).hexdigest()


def sha256_path(path: Path) -> str:
    h = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def safe_extract(zip_path: Path, destination: Path) -> None:
    with zipfile.ZipFile(zip_path) as zf:
        for info in zf.infolist():
            candidate = Path(info.filename)
            if candidate.is_absolute() or ".." in candidate.parts:
                raise RuntimeError(f"unsafe realm archive member: {info.filename!r}")
        zf.extractall(destination)


def normalize_label(value: object) -> str:
    return " ".join(str(value).strip().lower().replace("_", " ").split())


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--zip", type=Path, required=True)
    ap.add_argument("--binding", type=Path, required=True)
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()

    binding = json.loads(args.binding.read_text())
    if binding.get("schema") != "ttf_genetic_palearctic_realm_source_binding_v0.1":
        raise RuntimeError("unexpected realm source binding")
    if binding.get("status") != "FROZEN_BEFORE_FORMAL_REALM_CLASSIFICATION":
        raise RuntimeError("realm source binding was not frozen before classification")

    expected_blob = binding["bound_processed_copy"]["vendored_git_blob_sha1"]
    observed_blob = git_blob_sha1(args.zip)
    if observed_blob != expected_blob:
        raise RuntimeError(
            f"vendored realm ZIP Git blob drift: {observed_blob} != {expected_blob}"
        )

    with tempfile.TemporaryDirectory(prefix="ttf_wwf_realms_") as td:
        root = Path(td)
        safe_extract(args.zip, root)
        shps = sorted(root.rglob("*.shp"))
        if len(shps) != 1:
            raise RuntimeError(f"expected exactly one .shp, found {len(shps)}")
        reader = shapefile.Reader(str(shps[0]))
        fields = [field[0] for field in reader.fields[1:]]
        records = list(reader.iterShapeRecords())

        expected_total = int(binding["formal_selection"]["expected_total_realm_records"])
        if len(records) != expected_total:
            raise RuntimeError(
                f"unexpected realm record count {len(records)} != {expected_total}"
            )

        accepted = {
            normalize_label(value)
            for value in binding["formal_selection"]["realm_label_accepted"]
        }
        matches = []
        inspected = []
        for index, sr in enumerate(records):
            attrs = {
                fields[i]: sr.record[i]
                for i in range(len(fields))
            }
            values = {normalize_label(value) for value in attrs.values()}
            inspected.append({
                "index": index,
                "attributes": {str(k): str(v) for k, v in attrs.items()},
            })
            if values & accepted:
                geom = sr.shape.__geo_interface__
                matches.append((index, attrs, geom))

        expected_pa = int(
            binding["formal_selection"]["expected_palearctic_records"]
        )
        if len(matches) != expected_pa:
            raise RuntimeError(
                "could not identify exactly one Palearctic realm record; "
                f"matched={len(matches)}, fields={fields}, records={inspected}"
            )

        index, attrs, geom = matches[0]
        feature = {
            "type": "Feature",
            "properties": {
                **{str(k): str(v) for k, v in attrs.items()},
                "realmcode": "PA",
                "formal_realm_name": "Palearctic",
                "source_record_index": index,
            },
            "geometry": geom,
        }

    payload = {
        "type": "FeatureCollection",
        "features": [feature],
        "ttf_provenance": {
            "binding_schema": binding["schema"],
            "binding_sha256": sha256_path(args.binding),
            "source_zip_sha256": sha256_path(args.zip),
            "source_zip_git_blob_sha1": observed_blob,
            "source_repository_blob_sha1": binding["bound_processed_copy"][
                "git_blob_sha1"
            ],
            "original_record_count": expected_total,
            "palearctic_record_count": 1,
            "response_used": False,
        },
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    print(json.dumps({
        "status": "PASS_BOUND_PALEARCTIC_REALM_CONVERSION",
        "source_zip_git_blob_sha1": observed_blob,
        "record_count": expected_total,
        "palearctic_records": 1,
        "fields": fields,
    }, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
