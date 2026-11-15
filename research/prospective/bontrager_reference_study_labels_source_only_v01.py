"""Bontrager primary-study source-label inventory: SOURCE METADATA ONLY.

Uses the previously SHA/MD5-verified Dryad original as hosted by Borealis.
The 47-column CSV is tokenized, but biological response columns (including
'fitness') are never selected, inspected, converted, exported, or modeled.
Do not interpret different published dataset DOIs or matching author-year
tokens as independent original experiments.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
import math
import re
import zipfile
from collections import Counter, defaultdict
from pathlib import Path
from urllib.request import Request, urlopen

ORIGINAL_URL = "https://borealisdata.ca/api/access/datafile/150984?format=original"
ORIGINAL_MEMBER = "dryad/data/fitness_for_analysis.csv"
ORIGINAL_BYTES = 811204
ORIGINAL_SHA256 = "2995d0f165b27ccf2f38f3d0b286845c107762a0c2b02df2290e9a5df6b396f0"
ORIGINAL_MD5 = "2de98c1bb2cd2ddffcb090eb65a74d9f"
ALLOWED_GROUPS = frozenset({"annual", "herbaceous_perennial", "woody_perennial"})
COLUMNS = ("study", "functional_group", "source", "site", "fitness_type",
           "taxon", "lat_source", "long_source", "lat_site", "long_site")
EXPECTED = {"studies": 103, "species": 107, "cells": 711}
YEAR = re.compile(r"\b(?:18|19|20)\d{2}\b")


def binomial(name: str) -> str | None:
    name = " ".join(str(name).strip().replace("_", " ").split())
    return name if re.fullmatch(r"[A-Z][A-Za-z-]+ [a-z][A-Za-z-]+", name) else None


def lonlat(lon: str, lat: str) -> tuple[float, float] | None:
    try:
        x, y = float(lon), float(lat)
    except (ValueError, TypeError, OverflowError):
        return None
    return (x, y) if (math.isfinite(x) and math.isfinite(y)
                      and -180 <= x <= 180 and -90 <= y <= 90) else None


def canonical(value: str) -> str:
    return " ".join(re.findall(r"[a-z0-9]+", value.casefold()))


def source_study_inventory(original_csv: bytes, expect: dict | None = None) -> dict:
    reader = csv.reader(
        io.TextIOWrapper(io.BytesIO(original_csv), encoding="utf-8-sig", newline=""),
        delimiter=",", strict=True,
    )
    headers = next(reader, [])
    if len(headers) != 47 or len(set(headers)) != 47 or not set(COLUMNS).issubset(headers):
        raise ValueError("Unrecognized original source header: must match 47-field CSV")
    if "fitness" not in headers or "fitness" in COLUMNS:
        raise ValueError("Fitness response must exist but must never be allowlisted")
    at = {key: headers.index(key) for key in COLUMNS}
    rows = 0
    source_coordinates = defaultdict(set)
    garden_coordinates = defaultdict(set)
    candidates = []
    for record in reader:
        if not record:
            continue
        rows += 1
        if len(record) != len(headers):
            raise ValueError("Original CSV row width drift")
        m = {key: record[position].strip() for key, position in at.items()}
        if m["functional_group"] not in ALLOWED_GROUPS:
            continue
        study = m["study"]
        sp = binomial(m["taxon"])
        if not (study and sp and m["source"] and m["site"] and m["fitness_type"]):
            continue
        origin = lonlat(m["long_source"], m["lat_source"])
        garden = lonlat(m["long_site"], m["lat_site"])
        if origin is None or garden is None:
            continue
        source_coordinates[(study, sp, m["source"])].add(origin)
        garden_coordinates[(study, sp, m["site"])].add(garden)
        candidates.append((study, sp, m["source"], m["site"], m["fitness_type"], origin))
    bad_sources = {key for key, values in source_coordinates.items() if len(values) != 1}
    bad_gardens = {key for key, values in garden_coordinates.items() if len(values) != 1}
    sources_by_cell = defaultdict(set)
    coords_by_cell = defaultdict(set)
    for study, sp, source, site, performance_type, origin in candidates:
        if (study, sp, source) in bad_sources or (study, sp, site) in bad_gardens:
            continue
        key = (study, sp, site, performance_type)
        sources_by_cell[key].add(source)
        coords_by_cell[key].add(origin)
    qualified = sorted(
        key for key, source_ids in sources_by_cell.items()
        if len(source_ids) >= 2 and len(coords_by_cell[key]) >= 2
    )
    taxon_by_study = defaultdict(set)
    cell_by_study = Counter()
    taxa = set()
    for study, sp, site, performance_type in qualified:
        taxon_by_study[study].add(sp)
        cell_by_study[study] += 1
        taxa.add(sp)
    if expect is not None and (
        len(taxon_by_study) != expect["studies"]
        or len(taxa) != expect["species"]
        or len(qualified) != expect["cells"]
    ):
        raise ValueError("Original frozen 103-study/107-plant/711-cell source geometry not reproduced")
    normalized = defaultdict(set)
    catalog = []
    for study in sorted(taxon_by_study):
        norm = canonical(study)
        normalized[norm].add(study)
        candidate_years = YEAR.findall(study)
        catalog.append({
            "original_study_label": study,
            "normalized_study_label": norm,
            "original_study_label_sha256": hashlib.sha256(study.encode("utf-8")).hexdigest(),
            "candidate_year_from_label_NOT_verified": candidate_years[0] if len(candidate_years) == 1 else None,
            "primary_paper_doi": None,
            "primary_paper_title": None,
            "qualified_plant_species": sorted(taxon_by_study[study]),
            "qualified_same_garden_cells": cell_by_study[study],
            "independent_original_experiment_certified": False,
        })
    collisions = [sorted(group) for group in normalized.values() if len(group) > 1]
    if collisions:
        raise ValueError("Different original study labels canonicalize to same key")
    return {
        "schema": "ttf_bontrager_primary_study_reference_labels_SOURCE_ONLY_v01",
        "status": "REFERENCE_CATALOG_PREPARED_NOT_STUDY_OVERLAP_MEASURED",
        "source_dataset_doi": "10.5061/dryad.pnvx0k6mc",
        "source_mirror_dataset_doi": "10.5683/SP2/DU3ZTQ",
        "source_original_member": ORIGINAL_MEMBER,
        "original_csv_physical_rows": rows,
        "distinct_qualified_original_study_labels": len(taxon_by_study),
        "distinct_qualified_plant_species": len(taxa),
        "qualified_study_species_garden_fitness_type_cells": len(qualified),
        "labels_array_sha256": hashlib.sha256("\n".join(sorted(taxon_by_study)).encode()).hexdigest(),
        "unresolved_original_primary_paper_dois": len(catalog),
        "outcome_columns_selected": [],
        "fitness_rows_fitted": 0,
        "Hargreaves_overlap_checked": False,
        "certified_independent_additional_studies": 0,
        "catalog": catalog,
    }


def run() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--output", type=Path, required=True)
    args = p.parse_args()
    req = Request(ORIGINAL_URL, headers={"User-Agent": "TTF-original-study-label-metadata-only/2026-11-15"})
    with urlopen(req, timeout=90) as remote:
        raw = remote.read(ORIGINAL_BYTES + 1)
    if (len(raw) != ORIGINAL_BYTES
        or hashlib.sha256(raw).hexdigest() != ORIGINAL_SHA256
        or hashlib.md5(raw).hexdigest() != ORIGINAL_MD5):
        raise ValueError("Source archive identity does not match prior frozen original")
    with zipfile.ZipFile(io.BytesIO(raw)) as archive:
        if ORIGINAL_MEMBER not in archive.namelist():
            raise ValueError("Source CSV missing")
        with archive.open(ORIGINAL_MEMBER) as reader:
            csv_bytes = reader.read(25_000_001)
    if len(csv_bytes) > 25_000_000:
        raise ValueError("Original source CSV too large")
    result = source_study_inventory(csv_bytes, EXPECTED)
    result["original_zip_SHA256_verified"] = ORIGINAL_SHA256
    result["original_zip_MD5_verified"] = ORIGINAL_MD5
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2, ensure_ascii=False, sort_keys=True) + "\n")
    print("BONTRAGER_SOURCE_LABEL_REFERENCE", result["status"],
          result["distinct_qualified_original_study_labels"],
          result["distinct_qualified_plant_species"],
          result["qualified_study_species_garden_fitness_type_cells"])


if __name__ == "__main__":
    run()
