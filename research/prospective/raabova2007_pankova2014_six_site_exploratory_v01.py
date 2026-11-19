"""Separate exploratory same-code six-site geography check against Pánková 2014.

The 2007 Raabová K1–K3, S1–S3 codes are retained in the original Hargreaves
transplant source and destination IDs; Pánková et al. (2014) PLOS ONE Table 1
publishes those exact six locality codes and DMS. This is a later paper,
not independent original 2007 field GPS and not the erroneous 2010/2008 crosswalk.
No fitness, climate, gene or speciation response is opened or modeled.
"""
from __future__ import annotations
import argparse
import hashlib
import json
import math
import re
from pathlib import Path

from bontrager_hargreaves_finegrain_overlap_v01 import original_bytes, haversine
from hargreaves_source_as_garden_georecovery_v01 import source_metadata, one_point
from raabova_2007_2008_original_ID_metadata_v01 import classify_code

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
CONTRACT = HERE / "raabova2007_pankova2014_six_site_exploratory_contract_v01.json"


def strict_dms(value: str, axis: str) -> float:
    m = re.fullmatch(r"(\d{1,3}):(\d{2}):(\d{2}(?:\.\d+)?)([NSEW])", value)
    if not m:
        raise ValueError("Expected original publisher whole-second DMS")
    deg, minute, seconds, hem = m.groups()
    deg, minute, seconds = int(deg), int(minute), float(seconds)
    if (minute >= 60 or seconds >= 60
        or axis == "lat" and (hem not in ("N", "S") or deg > 90)
        or axis == "lon" and (hem not in ("E", "W") or deg > 180)
        or axis not in ("lon", "lat")):
        raise ValueError("Invalid original published DMS")
    v = deg + minute/60 + seconds/3600
    return -v if hem in ("W", "S") else v


def fixed_published_sites(c: dict) -> dict:
    primary = c["target"]
    other = c["external_primary_author_source"]
    codes = primary["source_population_codes"]
    if codes != ["K1", "K2", "K3", "S1", "S2", "S3"]:
        raise ValueError("Frozen original six population codes changed")
    if other["doi"] != "10.1371/journal.pone.0093967":
        raise ValueError("Not the original PLOS published site coordinates")
    places = other["published_site_rows"]
    if [p["code"] for p in places] != codes or len(places) != 6:
        raise ValueError("Frozen six published site rows changed")
    for p in places:
        if p["region"] != p["code"][0]:
            raise ValueError("Study-region and site code disagree")
    return {p["code"]: (
        strict_dms(p["lon_dms"], "lon"), strict_dms(p["lat_dms"], "lat")
    ) for p in places}


def coded_garden_proxy_coordinates(metadata: dict, c: dict) -> dict:
    primary = c["target"]
    accepted = primary["source_population_codes"]
    study, taxon = primary["original_study"], primary["strict_species"]
    observed = [key for key in metadata["origin_keys"]
                if key[0] == study and key[1] == taxon]
    if len(observed) != primary["expected_original_unique_source_codes"]:
        raise ValueError("Original 2007 source population count drift")
    by_code = {}
    for key in observed:
        code = classify_code(key[2], accepted)
        if code not in accepted or code in by_code:
            raise ValueError("Ambiguous or nonunique original K/S source ID")
        garden = metadata["gardens"].get(key, set())
        coords = one_point(garden)
        if coords is None:
            raise ValueError("No unique original Hargreaves destination coordinate for code")
        by_code[code] = coords
    if set(by_code) != set(accepted):
        raise ValueError("Original K/S population code completeness lost")
    return by_code


def exploratory_comparison(metadata: dict, contract: dict, erratum: dict) -> dict:
    if contract["status"] != "FROZEN_BEFORE_HARGREAVES_CODE_MATCHED_SIX_SITE_GEO_COMPARISON":
        raise ValueError("Six-site source-only design was not frozen")
    if (erratum["status"] != "TERMINAL_INVALID_EXTERNAL_VALIDATION_WRONG_PRIMARY_PAPER_REFERENT"
        or erratum["author_publication_referent"]["referenced_primary_DOI"] != "10.1007/s00442-008-1156-1"):
        raise ValueError("Previous invalid referent must remain corrected")
    r = contract["rules"]
    if (r["precommitted_descriptive_distance_bins_km"] != [0.1, 1, 2, 5, 10]
        or not r["match_six_original_population_codes_by_literal_code_not_geographically_nearest"]
        or not r["descriptive_previous_1km_line_not_a_pass_gate_for_new_exploration"]
        or not contract["no_biological_response_opening"]):
        raise ValueError("Source-only exploratory matching/reporting rules drifted")
    pub = fixed_published_sites(contract)
    proxy = coded_garden_proxy_coordinates(metadata, contract)
    results = []
    for code in contract["target"]["source_population_codes"]:
        distance_km = haversine(pub[code], proxy[code])
        cutoffs = r["precommitted_descriptive_distance_bins_km"]
        if distance_km <= cutoffs[0]:
            distance_bin = "0-0.1km"
        elif distance_km <= cutoffs[1]:
            distance_bin = ">0.1-1km"
        elif distance_km <= cutoffs[2]:
            distance_bin = ">1-2km"
        elif distance_km <= cutoffs[3]:
            distance_bin = ">2-5km"
        elif distance_km <= cutoffs[4]:
            distance_bin = ">5-10km"
        else:
            distance_bin = ">10km"
        results.append({
            "original_2007_population_code":code,
            "published_later_population_name":next(
                p["name"] for p in contract["external_primary_author_source"]["published_site_rows"]
                if p["code"] == code),
            "region":code[0],
            "distance_km":round(distance_km, 6),
            "distance_bucket_predeclared":distance_bin,
            "descriptively_le_1km":distance_km <= 1,
            "original_Hargreaves_source_population_ID_or_coordinates_exported":False
        })
    km = [x["distance_km"] for x in results]
    return {
        "schema":"ttf_Raabova2007_Pankova2014_same_code_six_site_exploratory_source_only_v01",
        "status":"NEW_EXPLORATORY_PUBLISHED_SIX_SAME_CODE_SITE_DISTANCES_NOT_FIELD_GPS_ACCURACY",
        "target_primary_paper_doi":contract["target"]["original_paper_doi"],
        "external_later_paper_doi":contract["external_primary_author_source"]["doi"],
        "exact_code_site_count":len(results),
        "same_region_site_codes_all_six_consistent":True,
        "site_details":results,
        "descriptive_max_distance_km":max(km),
        "descriptive_sites_within_1km":sum(x["descriptively_le_1km"] for x in results),
        "descriptive_sites_over_1km":sum(not x["descriptively_le_1km"] for x in results),
        "source_garden_proxy_recovered_by_within_study_ID_not_spatial_nearest":True,
        "external_paper_is_2014_followup_NOT_independent_field_GPS":True,
        "old_invalid_Sudova_2010_vs_2007_receipt_remains_invalid":True,
        "no_original_source_ID_or_raw_dataset_point_coordinates_exported":True,
        "independent_population_GPS_validation_certified":False,
        "independent_primary_experiment_certified":False,
        "no_adaptation_plasticity_speciation_genetic_history_climate_response_opened":True,
        "old_Bontrager_44_of_50_source_gate_not_requalified":True
    }


def run():
    p = argparse.ArgumentParser()
    p.add_argument("--output", required=True, type=Path)
    args = p.parse_args()
    c = json.loads(CONTRACT.read_text(encoding="utf-8"))
    errs = json.loads((ROOT/c["hargreaves_original"]["correct_2007_vs_2008_erratum"]).read_text())
    orig = json.loads((ROOT/c["hargreaves_original"]["original_source_contract"]).read_text())
    entry = orig["original"]["hargreaves"]
    raw = original_bytes(entry)
    if hashlib.sha256(raw).hexdigest() != c["hargreaves_original"]["csv_sha256"]:
        raise ValueError("Original Hargreaves input identity changed")
    r = json.loads((HERE/"hargreaves_source_as_garden_georecovery_contract_v01.json").read_text())
    meta = source_metadata(raw, "hargreaves", r)
    out = exploratory_comparison(meta, c, errs)
    out["original_Hargreaves_CSV_SHA256_verified"] = c["hargreaves_original"]["csv_sha256"]
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(out, indent=2, sort_keys=True, ensure_ascii=False)+"\n", encoding="utf-8")
    print("RAABOVA_PANKOVA_SIX_SITE_SOURCE_QA=", json.dumps({
        "codes": [x["original_2007_population_code"] for x in out["site_details"]],
        "km": [x["distance_km"] for x in out["site_details"]],
        "le_1km":out["descriptive_sites_within_1km"],
        "status":out["status"]}, sort_keys=True))

if __name__ == "__main__":
    run()
