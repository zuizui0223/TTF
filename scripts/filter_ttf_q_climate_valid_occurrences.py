#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from collections import Counter
from pathlib import Path

import numpy as np

from ttf.climate_occurrence_filter import climate_validity_mask


def sha256_path(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--occurrences", type=Path, required=True)
    ap.add_argument("--raster", type=Path, action="append", required=True)
    ap.add_argument("--output", type=Path, required=True)
    ap.add_argument("--summary", type=Path, required=True)
    args = ap.parse_args()

    if len(args.raster) < 1:
        raise RuntimeError("at least one climate raster is required")

    rows = list(csv.DictReader(args.occurrences.open(encoding="utf-8")))
    if not rows:
        raise RuntimeError("empty occurrence table")
    coords = [
        (float(row["longitude"]), float(row["latitude"]))
        for row in rows
    ]

    import rasterio

    datasets = [rasterio.open(path) for path in args.raster]
    try:
        columns = []
        nodata = []
        for dataset in datasets:
            if dataset.crs is None or not dataset.crs.is_geographic:
                raise RuntimeError(
                    f"climate raster is not geographic: {dataset.name}"
                )
            values = np.asarray(
                [float(sample[0]) for sample in dataset.sample(coords)],
                dtype=float,
            )
            columns.append(values)
            nodata.append(
                None if dataset.nodata is None else float(dataset.nodata)
            )
        matrix = np.column_stack(columns)
        valid = climate_validity_mask(matrix, nodata_values=nodata)
    finally:
        for dataset in datasets:
            dataset.close()

    retained = [row for row, keep in zip(rows, valid) if bool(keep)]
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(retained)

    before = Counter(row["species"] for row in rows)
    after = Counter(row["species"] for row in retained)
    payload = {
        "schema": "ttf_q_climate_valid_occurrence_filter_v0.1",
        "status": "FILTERED_RESPONSE_BLIND_CLIMATE_VALID_OCCURRENCES",
        "input_rows": len(rows),
        "retained_rows": len(retained),
        "removed_rows": len(rows) - len(retained),
        "input_species": len(before),
        "retained_species": len(after),
        "species_ge_30_before": sum(n >= 30 for n in before.values()),
        "species_ge_30_after": sum(n >= 30 for n in after.values()),
        "input_occurrence_sha256": sha256_path(args.occurrences),
        "output_occurrence_sha256": sha256_path(args.output),
        "rasters": {
            path.name: sha256_path(path)
            for path in args.raster
        },
        "genetic_response_used": False,
        "purpose": (
            "Prevent finite/nodata-invalid raster cells from entering "
            "historical displacement arithmetic in fresh TTF-Q method development."
        ),
    }
    args.summary.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n"
    )
    print(json.dumps({
        "status": payload["status"],
        "retained_rows": payload["retained_rows"],
        "species_ge_30_after": payload["species_ge_30_after"],
    }, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
