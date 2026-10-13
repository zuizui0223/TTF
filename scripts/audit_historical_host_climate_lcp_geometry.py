#!/usr/bin/env python3
"""Describe already observed CHELSA LCP edge costs without using any genetic response.

The original 26-locality, 64-edge exploratory CHELSA LCP result is immutable.
This descriptive audit examines scale and sampling geometry; it is not a new
hypothesis test, ancient plant-distribution model or evidence of migration.
"""
from __future__ import annotations
import argparse
import csv
import hashlib
import json
import math
from pathlib import Path

import numpy as np

SPECIES = "Bungalotis erythus"
ORIGINAL_SUMMARY_SHA256 = "6969f482a1967db2db159cbc5fd62b722268a0fc45d7dcd1b979f1973b9a3493"
ORIGINAL_EDGE_SHA256 = "7d255b208a485e7bc6d02e4e5748781c5e10170ed8bf161542dce3135c897b3b"
SOURCE_ARTIFACT_ZIP_SHA256 = "0a955188675f98e1597e6b2f59b7ba31e09aa6b640f745a0459976cb06c50aa3"


def sha(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for b in iter(lambda:f.read(1 << 20),b""):
            h.update(b)
    return h.hexdigest()


def audit_rows(rows: list[dict], *, expected_species: str, expected_edges: int) -> dict:
    if len(rows) != expected_edges or not rows:
        raise ValueError("frozen edge count drift")
    if {row.get("species") for row in rows} != {expected_species}:
        raise ValueError("frozen source species mismatch")
    ids = [int(row["edge_index"]) for row in rows]
    if sorted(ids) != list(range(expected_edges)):
        raise ValueError("original frozen edge indices changed")
    edges = [(int(row["node_left"]),int(row["node_right"])) for row in rows]
    if any(l < 0 or r < 0 or l == r for l,r in edges):
        raise ValueError("invalid original graph edge")
    unique = [tuple(sorted(pair)) for pair in edges]
    if len(set(unique)) != len(unique):
        raise ValueError("duplicate undirected frozen edge")
    keys = ("geographic_chord_km","R_host_0BP_km","R_host_LGM_km",
            "R_self_0BP_km","R_self_LGM_km","M_host_lcp_km",
            "M_self_lcp_km","M_host_minus_self_lcp_km")
    vals = {k:np.asarray([float(row[k]) for row in rows],dtype=float) for k in keys}
    if any(not np.all(np.isfinite(a)) for a in vals.values()):
        raise ValueError("nonfinite edge cost / original geography")
    chord = vals["geographic_chord_km"]
    if np.any(chord <= 0):
        raise ValueError("original chord distance must be positive")
    for k in ("R_host_0BP_km","R_host_LGM_km","R_self_0BP_km","R_self_LGM_km"):
        if np.any(vals[k] < 0):
            raise ValueError("negative path costs")
    dh = vals["M_host_lcp_km"]
    ds = vals["M_self_lcp_km"]
    dd = vals["M_host_minus_self_lcp_km"]
    if (np.max(np.abs(dh - (vals["R_host_LGM_km"] - vals["R_host_0BP_km"]))) > 1e-11 or
        np.max(np.abs(ds - (vals["R_self_LGM_km"] - vals["R_self_0BP_km"]))) > 1e-11 or
        np.max(np.abs(dd - (dh - ds))) > 1e-11):
        raise ValueError("LCP edge delta does not match four frozen scenarios")
    same_pixel = ((vals["R_host_0BP_km"] == 0) & (vals["R_host_LGM_km"] == 0) &
                  (vals["R_self_0BP_km"] == 0) & (vals["R_self_LGM_km"] == 0))
    if np.std(chord) <= 1e-12 or np.std(dh) <= 1e-12 or np.std(ds) <= 1e-12:
        raise ValueError("cannot estimate descriptive chord correlation for invariant test input")
    pearson = lambda a,b:float(np.corrcoef(a,b)[0,1])
    return {
        "schema":"ttf_historical_host_climate_native_LCP_geometry_diagnostic_v0.1",
        "status":"EXPLORATORY_CLIMATE_ONLY_GEO_NORMALIZATION",
        "species":expected_species,
        "n_original_edges":int(len(rows)),
        "counts":{
            "host_LGM_cost_higher":int(np.sum(dh>0)),
            "host_LGM_cost_equal":int(np.sum(dh==0)),
            "self_LGM_cost_higher":int(np.sum(ds>0)),
            "self_LGM_cost_equal":int(np.sum(ds==0)),
            "host_minus_self_negative":int(np.sum(dd<0)),
            "host_minus_self_equal":int(np.sum(dd==0)),
            "same_cell_zero_path_edges":int(np.sum(same_pixel)),
        },
        "pairwise_km_normalized":{
            "host_mean_delta_per_original_chord_km":float(np.mean(dh/chord)),
            "self_mean_delta_per_original_chord_km":float(np.mean(ds/chord)),
            "host_minus_self_mean_delta_per_original_chord_km":float(np.mean(dd/chord)),
            "host_median_delta_per_original_chord_km":float(np.median(dh/chord)),
            "self_median_delta_per_original_chord_km":float(np.median(ds/chord)),
        },
        "geometry_association":{
            "pearson_original_chord_vs_host_delta":pearson(chord,dh),
            "pearson_original_chord_vs_self_delta":pearson(chord,ds),
            "pearson_original_chord_vs_host_minus_self_delta":pearson(chord,dd),
        },
        "diagnostics":{
            "no_genetic_effect_or_p_value":True,
            "no_paleo_land_ice_mask":True,
            "no_new_genetic_response_read":True,
            "all_64_edges_share_26_locations_not_independent":True,
            "cloud_support_size_asymmetry_makes_cross_reference_resistance_not_directly_comparable":True,
            "original_climate_cost_metric_not_retuned":True,
            "one_same_raster_cell_path_zero_is_resolution_artifact":True,
            "distance_normalization_is_descriptive_not_an_alternative_hypothesis_test":True,
            "modern_host_niche_is_not_proven_LGM_plant_occupancy":True,
            "original_259_confirmatory_status":"NOT_EVALUABLE_EMPIRICAL_RESPONSE",
        },
    }


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--pilot-summary",required=True,type=Path)
    p.add_argument("--pilot-edges",required=True,type=Path)
    p.add_argument("--output",required=True,type=Path)
    a = p.parse_args()
    if sha(a.pilot_summary) != ORIGINAL_SUMMARY_SHA256 or sha(a.pilot_edges) != ORIGINAL_EDGE_SHA256:
        raise RuntimeError("original immutable 64-edge CHELSA pilot source SHA mismatch")
    summary = json.loads(a.pilot_summary.read_text())
    if (summary["status"] != "PASS_CHELSA_CLIMATE_ONLY_PILOT_NO_LAND_ICE_MASK" or
        summary["species"] != SPECIES or summary["localities"] != 26 or
        summary["original_edges"] != 64 or
        summary["output_edges_sha256"] != ORIGINAL_EDGE_SHA256 or
        max(summary["frozen_5point_predictor_max_abs_error"].values()) != 0.0):
        raise RuntimeError("frozen original climate-only pilot status or 5-point replay drift")
    with a.pilot_edges.open(newline="", encoding="utf-8") as f:
        rows = list(csv.DictReader(f))
    result = audit_rows(rows,expected_species=SPECIES,expected_edges=64)
    if result["counts"]["same_cell_zero_path_edges"] != summary["same_cell_original_edges"]:
        raise RuntimeError("same pixel path ledger drift")
    if result["counts"]["host_LGM_cost_higher"] != 63 or result["counts"]["host_minus_self_negative"] != 63:
        raise RuntimeError("original source 64-edge signed counts drift")
    result.update({
        "input_sha256": {"pilot_summary":ORIGINAL_SUMMARY_SHA256,"pilot_edges":ORIGINAL_EDGE_SHA256,"original_artifact_zip":SOURCE_ARTIFACT_ZIP_SHA256},
        "host_0bp_climate_cloud_points": summary["host_0bp_climate_cloud_points"],
        "self_0bp_climate_cloud_points": summary["self_0bp_climate_cloud_points"],
        "same_cell_original_edge_ids": [int(row["edge_index"]) for row in rows if (float(row["R_host_0BP_km"]) == 0 and float(row["R_self_LGM_km"]) == 0)],
        "original_aggregate_mean_delta_host": summary["change_LGM_minus_present_cost_km"]["host"]["mean"],
        "original_aggregate_mean_delta_self": summary["change_LGM_minus_present_cost_km"]["insect_self"]["mean"],
    })
    if result["host_0bp_climate_cloud_points"] != 665 or result["self_0bp_climate_cloud_points"] != 26:
        raise RuntimeError("fixed host/self climate cloud reference support sizes drift")
    a.output.parent.mkdir(parents=True,exist_ok=True)
    a.output.write_text(json.dumps(result,sort_keys=True,indent=2)+"\n",encoding="utf-8")
    print(json.dumps({"status":result["status"],"n":result["n_original_edges"],"pearson":result["geometry_association"],"per_km":result["pairwise_km_normalized"]},sort_keys=True))


if __name__ == "__main__":
    main()
