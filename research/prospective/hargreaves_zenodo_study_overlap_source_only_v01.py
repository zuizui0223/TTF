"""Response-blind Hargreaves 2020 original source study-label overlap against Bontrager.

Only source study, taxon, kingdom, functional group, source and garden labels
are selected from the 47 original CSV fields. Original full bytes are hashed
against publisher Zenodo MD5 before ANY row is parsed. Original biological
fitness, treatment effects, outcomes and location-derived response fields
are NOT indexed, converted to numbers, examined or modeled.

Exact original source-study labels are strong provenance overlap candidates,
not proof that two different meta-analysis rows refer to the same experiment.
Different/unmatched labels are NOT certified independent experiments.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
import re
from collections import Counter, defaultdict
from pathlib import Path
from urllib.request import Request, urlopen

from hargreaves_zenodo_header_source_only_v01 import (
    RECORD_URL, original_file_info
)

HERE = Path(__file__).resolve().parent
CONTRACT_PATH = HERE / "hargreaves_zenodo_study_overlap_contract_v01.json"
MAX_ROWS = 200_000


def canonical(value: str) -> str:
    return " ".join(re.findall(r"[a-z0-9]+", value.casefold().replace("_", " ")))


def source_only_study_overlap(original_csv: bytes, contract: dict, reference: dict) -> dict:
    source = contract["source"]
    if (len(original_csv) != source["public_file_metadata_bytes"]
        or hashlib.md5(original_csv).hexdigest() != source["public_file_metadata_md5"]):
        raise ValueError("Original Zenodo CSV bytes/md5 mismatch; no source rows authorized")
    if b"\n" not in original_csv:
        raise ValueError("Original CSV header line missing")
    header_bytes = original_csv.split(b"\n", 1)[0] + b"\n"
    if hashlib.sha256(header_bytes).hexdigest() != source["verified_header_firstline_sha256"]:
        raise ValueError("Original source-only header SHA changed")
    # Source encoding is decided by byte-level validity BEFORE any source ID
    # or taxon value is indexed. The original data file is never rewritten.
    try:
        decoded = original_csv.decode("utf-8-sig", errors="strict")
        selected_encoding = "utf-8-sig"
    except UnicodeDecodeError:
        decoded = original_csv.decode("cp1252", errors="strict")
        selected_encoding = "cp1252"
    reader = csv.reader(io.StringIO(decoded, newline=""), strict=True)
    header = next(reader)
    if header != source["original_csv_header"] or len(header) != 47:
        raise ValueError("Source original 47-column schema not identical to precommitted header")
    allowed = contract["comparability"]["metadata_columns_allowlist"]
    if (allowed != ["study", "functional_group", "kingdom", "taxon", "source", "site"]
        or any(x in allowed for x in ("fitness","fit_var","thermal_dist","lat_site","long_site"))):
        raise ValueError("Source-only metadata allowlist drifted into responses or source coordinates")
    positions = {name: header.index(name) for name in allowed}

    base = contract["base_reference"]
    if (len(reference["catalog"]) != base["original_study_label_count"]
        or reference["labels_array_sha256"] != "83a80a48be293a3ca2b91ebafcea6e4860ba5c7100c7d11447e3aa71c107fa85"
        or reference["distinct_qualified_plant_species"] != base["original_plant_species_count"]
        or reference["qualified_study_species_garden_fitness_type_cells"] != base["original_shared_garden_cells"]):
        raise ValueError("Original source reference geometry drifted")

    plant_groups = set(contract["comparability"]["original_plant_functional_groups"])
    all_studies, plants, plant_studies = set(), set(), set()
    taxa_per_study = defaultdict(set)
    functional_group_counts, kingdom_counts = Counter(), Counter()
    total_rows, plant_rows = 0, 0
    for row in reader:
        if not row:
            continue
        total_rows += 1
        if total_rows > MAX_ROWS or len(row) != len(header):
            raise ValueError("Original source CSV row count or width outside bounded contract")
        m = {key: row[idx].strip() for key, idx in positions.items()}
        functional_group_counts[m["functional_group"] or "<empty>"] += 1
        kingdom_counts[m["kingdom"] or "<empty>"] += 1
        study = m["study"]
        if study:
            all_studies.add(study)
        if m["functional_group"] not in plant_groups:
            continue
        plant_rows += 1
        if study:
            plant_studies.add(study)
        taxon = canonical(m["taxon"])
        if taxon:
            plants.add(taxon)
            if study:
                taxa_per_study[study].add(taxon)

    ref_by_label = {
        item["original_study_label"]: {
            "taxa": {canonical(sp) for sp in item["qualified_plant_species"]},
        } for item in reference["catalog"]
    }
    norm_to_ref = defaultdict(set)
    for label in ref_by_label:
        norm_to_ref[canonical(label)].add(label)

    exact, normalized_candidates, unmatched = [], [], []
    exact_plant = []
    shared_study_taxon_pairs = 0
    for label in sorted(plant_studies):
        if label in ref_by_label:
            exact.append(label)
            same_taxa = taxa_per_study[label] & ref_by_label[label]["taxa"]
            shared_study_taxon_pairs += len(same_taxa)
            if same_taxa:
                exact_plant.append(label)
        else:
            aliases = norm_to_ref.get(canonical(label), set())
            if len(aliases) == 1:
                normalized_candidates.append(label)
            else:
                unmatched.append(label)
    return {
        "schema": "ttf_Hargreaves2020_vs_Bontrager_original_study_SOURCE_ONLY_overlap_v01",
        "status": "ORIGINAL_SOURCE_STUDY_ID_OVERLAP_MEASURED_NOT_INDEPENDENT_EXPERIMENT_CERTIFICATION",
        "source_original_CSV_bytes": len(original_csv),
        "source_original_CSV_md5": hashlib.md5(original_csv).hexdigest(),
        "source_original_CSV_sha256": hashlib.sha256(original_csv).hexdigest(),
        "source_original_text_decode_selected_after_byte_md5_check": selected_encoding,
        "raw_original_bytes_not_reencoded_or_modified": True,
        "original_source_physical_rows": total_rows,
        "original_plant_functional_group_rows": plant_rows,
        "original_source_functional_group_counts": dict(sorted(functional_group_counts.items())),
        "original_source_kingdom_counts": dict(sorted(kingdom_counts.items())),
        "hargreaves_all_source_study_labels": len(all_studies),
        "hargreaves_plant_source_study_labels": len(plant_studies),
        "hargreaves_plant_source_taxa_text_labels": len(plants),
        "Bontrager_reference_original_study_labels": len(ref_by_label),
        "exact_original_study_label_overlap_plant": len(exact),
        "exact_source_label_and_plant_taxon_joint_keys": shared_study_taxon_pairs,
        "exact_original_study_labels_with_at_least_one_same_taxon": len(exact_plant),
        "normalized_only_source_study_label_candidates_for_review": len(normalized_candidates),
        "unlinked_source_labels_NOT_independence_certified": len(unmatched),
        "exact_source_study_label_overlap_labels": sorted(exact),
        "normalized_only_source_study_labels_review": sorted(normalized_candidates),
        "unlinked_source_study_labels_not_certified": sorted(unmatched),
        "source_ID_and_garden_ID_values_exported": False,
        "source_original_response_and_effect_columns_used": [],
        "source_original_numeric_fitness_rows_fit": 0,
        "full_original_bytes_hashed_before_source_label_parsing": True,
        "no_response_or_adaptive_fitness_model_opened": True,
        "source_overlap_NOT_experiment_identity_certification": True,
        "unlinked_not_proof_of_independent_primary_experiment": True,
        "Bontrager_frozen_44_of_50_gate_unchanged": True,
    }


def run() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--output", required=True, type=Path)
    args = ap.parse_args()
    contract = json.loads(CONTRACT_PATH.read_text(encoding="utf-8"))
    if contract["status"] != "FROZEN_BEFORE_ORIGINAL_HARGREAVES_DATA_ROW_METADATA_CENSUS":
        raise ValueError("No pre-response original metadata-only contract")
    base = contract["base_reference"]
    raw_ref = (HERE.parents[1] / base["verified_reference_json_path"]).read_bytes()
    if hashlib.sha256(raw_ref).hexdigest() != base["reference_json_sha256"]:
        raise ValueError("Source reference Bontrager catalog SHA drifted")
    reference = json.loads(raw_ref)
    with urlopen(Request(RECORD_URL,headers={"Accept":"application/json",
       "User-Agent":"TTF-preserved-Hargreaves-source-only-study-audit/2026-11-15"}),timeout=30) as remote:
        record_raw = remote.read(250_001)
    if len(record_raw) > 250_000:
        raise ValueError("Official Zenodo record metadata exceeds cap")
    info = original_file_info(json.loads(record_raw))
    md5_expected = "md5:" + contract["source"]["public_file_metadata_md5"]
    if (info["size"] != contract["source"]["public_file_metadata_bytes"]
        or (info["checksum_from_published_metadata"] or "").lower() != md5_expected):
        raise ValueError("Official Zenodo file identity metadata changed")
    with urlopen(Request(info["content_url"],headers={"Accept":"*/*",
        "User-Agent":"TTF-Hargreaves-source-label-census-only/2026-11-15"}),timeout=75) as remote:
        raw = remote.read(contract["source"]["public_file_metadata_bytes"] + 1)
    result = source_only_study_overlap(raw, contract, reference)
    result["official_published_zenodo_record_id"] = 4409291
    result["original_Dryad_DOI"] = contract["source"]["original_dryad_doi"]
    result["reference_source_label_catalog_SHA256"] = base["reference_json_sha256"]
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(result,indent=2,ensure_ascii=False,sort_keys=True)+"\n",encoding="utf-8")
    print("HARGREAVES_BONTRAGER_SOURCE_OVERLAP=", json.dumps({
        k: result[k] for k in (
            "status", "original_source_physical_rows", "hargreaves_all_source_study_labels",
            "hargreaves_plant_source_study_labels", "hargreaves_plant_source_taxa_text_labels",
            "exact_original_study_label_overlap_plant",
            "exact_source_label_and_plant_taxon_joint_keys",
            "normalized_only_source_study_label_candidates_for_review",
            "unlinked_source_labels_NOT_independence_certified",
        )
    },sort_keys=True))


if __name__ == "__main__":
    run()
