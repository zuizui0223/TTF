"""No genetics; frozen CHELSA native-grid pilot receipt and 64-edge verification."""
from __future__ import annotations
import csv
import hashlib
import json
import math
from pathlib import Path
import statistics

ROOT = Path(__file__).resolve().parents[1]
DIR = ROOT / "benchmarks" / "exploratory"
SOURCE_ZIP_SHA256 = "0a955188675f98e1597e6b2f59b7ba31e09aa6b640f745a0459976cb06c50aa3"
SUMMARY_SHA256 = "6969f482a1967db2db159cbc5fd62b722268a0fc45d7dcd1b979f1973b9a3493"
EDGES_SHA256 = "7d255b208a485e7bc6d02e4e5748781c5e10170ed8bf161542dce3135c897b3b"
SUMMARY = DIR / "historical_host_climate_lcp_bungalotis_summary_v0.1.json"
EDGES = DIR / "historical_host_climate_lcp_bungalotis_64edges_v0.1.csv"


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def test_preserve_exact_GHA_artifact_data_and_response_firewall():
    assert sha256(SUMMARY) == SUMMARY_SHA256
    assert sha256(EDGES) == EDGES_SHA256
    p = json.loads(SUMMARY.read_text())
    assert p["status"] == "PASS_CHELSA_CLIMATE_ONLY_PILOT_NO_LAND_ICE_MASK"
    assert p["species"] == "Bungalotis erythus"
    assert p["localities"] == 26 and p["original_edges"] == 64
    assert p["output_edges_sha256"] == EDGES_SHA256
    assert p["contract_sha256"] == "a46e8c61bbc07e3b093c4664b12825c7b193455be99d457d8d96343a358d5ef7"
    assert p["ecological_boundary"]["COI_genetic_response_used"] is False
    assert p["ecological_boundary"]["no_paleo_land_or_ice_mask_accepted"] is True
    assert p["ecological_boundary"]["no_genetic_effect_p_value"] is True
    assert max(p["frozen_5point_predictor_max_abs_error"].values()) == 0.0


def test_all_original_64_edges_and_paired_delta_arithmetic():
    s = json.loads(SUMMARY.read_text())
    with EDGES.open(newline="", encoding="utf-8") as f:
        rows = list(csv.DictReader(f))
    assert len(rows) == 64
    assert [int(row["edge_index"]) for row in rows] == list(range(64))
    keys = ("R_host_0BP_km", "R_host_LGM_km", "R_self_0BP_km", "R_self_LGM_km", "M_host_lcp_km", "M_self_lcp_km", "M_host_minus_self_lcp_km")
    for row in rows:
        vals = {k:float(row[k]) for k in keys}
        assert all(math.isfinite(v) and v >= 0 for k,v in vals.items() if k.startswith("R_"))
        assert math.isclose(vals["M_host_lcp_km"],vals["R_host_LGM_km"] - vals["R_host_0BP_km"],abs_tol=1e-12)
        assert math.isclose(vals["M_self_lcp_km"],vals["R_self_LGM_km"] - vals["R_self_0BP_km"],abs_tol=1e-12)
        assert math.isclose(vals["M_host_minus_self_lcp_km"],vals["M_host_lcp_km"] - vals["M_self_lcp_km"],abs_tol=1e-12)
    for kind,column in (("host","M_host_lcp_km"),("insect_self","M_self_lcp_km"),("host_minus_insect_self","M_host_minus_self_lcp_km")):
        x = [float(row[column]) for row in rows]
        z = s["change_LGM_minus_present_cost_km"][kind]
        assert math.isclose(sum(x)/len(x),z["mean"],abs_tol=1e-12)
        assert math.isclose(statistics.median(x),z["median"],abs_tol=1e-12)
        assert min(x) == z["min"] and max(x) == z["max"]
    assert sum(float(r["M_host_lcp_km"])>0 for r in rows) == 63
    assert sum(float(r["M_host_minus_self_lcp_km"])<0 for r in rows) == 63
    assert s["same_cell_original_edges"] == 1


def test_no_claim_of_ecological_causation_from_reference_cloud_size_or_raw_edge_count():
    s = json.loads(SUMMARY.read_text())
    assert s["host_0bp_climate_cloud_points"] == 665
    assert s["self_0bp_climate_cloud_points"] == 26
    assert s["native_raster"]["rows"] * s["native_raster"]["cols"] == s["joint_valid_cells"] == 5757
    assert s["ecological_boundary"]["modern_host_climate_analogue_is_not_ancient_host_occupancy"] is True
    assert s["ecological_boundary"]["ROI_boundary_can_artificially_limit_detours"] is True
    assert s["ecological_boundary"]["pilot_not_confirmatory_or_source_disjoint"] is True
    assert s["ecological_boundary"]["original_259_empirical_result"] == "NOT_EVALUABLE_EMPIRICAL_RESPONSE"
