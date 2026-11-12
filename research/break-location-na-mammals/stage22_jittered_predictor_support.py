"""Stage22: re-sample all 13 frozen predictors at 50-km-jittered Jensen cases.

Response-blind source/predictor audit ONLY. The 80-species original case sites,
surveyed backgrounds, sampling density, sources and SHA locks remain frozen.
No real genetic break model or AUC is read, fitted or calculated.
All IUCN-derived geometry and location matrices remain in ephemeral runner RAM.
"""
from __future__ import annotations

import argparse
import gc
import hashlib
import json
import platform
import zipfile
from pathlib import Path

import numpy as np
import pyproj
import scipy
import shapely
import rasterio
from pyproj import CRS, Transformer
from shapely.ops import transform
from shapely.validation import make_valid

from extract_frozen_site_features import (
    rebuild_sites, obtain, load, vector_geometries, mountain_distances,
)
from jitter_case_geometry import jitter_cases, repair_survey_polygon
from source_overlap_gate import members, load_group, norm
from site_feature_core import (
    complete_rows, dense_vertices, raster_values, vertex_distances,
)

# Authoritative original Stage20 sampling-geometry/case hashes. Never change.
ORIGINAL_XY_SHA256 = "fad49f6ec345ef45c1387ef76c8b3b724caf5f39322f8e9c3da56227a7a40bfd"
ORIGINAL_CASE_SHA256 = "e08632a14a990804a7a1c83e236b70834056ac61b7157b9832b8f91c0df855d5"
REPLAY_FEATURE_SHA256 = "51cf2953dda21ebf06f5cb78bbc0cf2c46f7a0d2bce1f7606128b27239bec2b8"


def sha_bytes(a, dtype):
    return hashlib.sha256(np.asarray(a).astype(dtype).tobytes()).hexdigest()


def ensure_frozen_geometry(xy, cases, species_index, records, panel):
    if len(xy) != 76700 or int(cases.sum()) != 10013 or len(records) != 80:
        raise ValueError("Stage20 80-species design changed")
    observed_xy_sha256 = sha_bytes(xy, "<f8")
    if observed_xy_sha256 != ORIGINAL_XY_SHA256:
        observed_case_sha256 = hashlib.sha256(np.asarray(cases, dtype="u1").tobytes()).hexdigest()
        raise ValueError("Stage20 coordinate identity STOP; expected=" + ORIGINAL_XY_SHA256 +
            " observed=" + observed_xy_sha256 + " case_labels=" + observed_case_sha256 +
            " numpy=" + np.__version__ + " shapely=" + shapely.__version__ +
            " pyproj=" + pyproj.__version__ + " (no biological response opened)")
    if hashlib.sha256(np.asarray(cases, dtype="u1").tobytes()).hexdigest() != ORIGINAL_CASE_SHA256:
        raise ValueError("Original case-order labels differ from frozen Stage20")
    roles = [sum(int(s["split_role"]) == r for s in records) for r in (0, 1)]
    if roles != [44, 36]:
        raise ValueError("Species heldout split changed")
    return roles


def jitter_cases_from_verified_sources(xy, cases, species_index, records, work, jitter_contract):
    """Recreate Stage21 local AEQD cases exactly; change ONLY case coordinates."""
    rules = jitter_contract["case_jitter"]
    if jitter_contract["status"] != "LOCKED_PRE_EMPIRICAL_BREAK_AUC":
        raise ValueError("Stage21 jitter protocol not frozen")
    if int(rules["source_uncertainty_radius_km"]) != 50 or int(rules["max_attempts"]) != 20:
        raise ValueError("50-km/20-draw position rule changed")
    moved_xy = np.array(xy, copy=True)
    case_shifts = []
    per_species = []
    with zipfile.ZipFile(work / "breaks.zip") as bz, zipfile.ZipFile(work / "ranges.zip") as rz:
        bb, rr = members(bz, "break"), members(rz, "range")
        for i, record in enumerate(records):
            key = norm(record["common_name"])
            line, _ = load_group(bz, bb[key][1])
            polygon, _ = load_group(rz, rr[key][1])
            if not polygon.is_valid:
                polygon = make_valid(polygon)
            anchor = polygon.representative_point()
            crs = CRS.from_proj4(
                f"+proj=aeqd +lat_0={anchor.y} +lon_0={anchor.x} "
                "+datum=WGS84 +units=m +no_defs"
            )
            forward = Transformer.from_crs(4326, crs, always_xy=True)
            inverse = Transformer.from_crs(crs, 4326, always_xy=True)
            projected_line = transform(forward.transform, line)
            projected_polygon = transform(forward.transform, polygon)
            d = np.arange(5000.0, projected_line.length, 10000.0)
            pts = shapely.line_interpolate_point(projected_line, d)
            pts = pts[shapely.covers(projected_polygon, pts)]
            local_xy = shapely.get_coordinates(pts)
            idx = np.flatnonzero((species_index == i) & (cases == 1))
            if len(idx) != int(record["case_points"]) or len(idx) != len(local_xy):
                raise ValueError("Original species case count differs: " + record["scientific_name"])
            old_lon, old_lat = inverse.transform(local_xy[:, 0], local_xy[:, 1])
            rebuilt_old = np.column_stack([old_lon, old_lat])
            if not np.allclose(rebuilt_old, xy[idx], rtol=0, atol=1e-8):
                raise ValueError("Stage21 AEQD case locations do not match frozen Stage20 cases")
            seed_key = rules["seed_namespace"] + record["scientific_name"].lower()
            seed = int.from_bytes(hashlib.sha256(seed_key.encode()).digest()[:8], "big")
            source_polygon = repair_survey_polygon(projected_polygon)
            moved, diagnostic = jitter_cases(
                local_xy, source_polygon, np.random.default_rng(seed),
                radius_m=rules["source_uncertainty_radius_km"] * 1000,
                max_attempts=rules["max_attempts"]
            )
            shift_km = np.linalg.norm(moved - local_xy, axis=1) / 1000.0
            if not shapely.covers(source_polygon, shapely.points(moved)).all():
                raise ValueError("Moved case left surveyed-range polygon")
            if (shift_km > 50.0 + 1e-8).any():
                raise ValueError("Position jitter exceeded 50 km")
            lon, lat = inverse.transform(moved[:, 0], moved[:, 1])
            moved_xy[idx] = np.column_stack([lon, lat])
            case_shifts.append(shift_km)
            per_species.append({
                "species": record["scientific_name"],
                "role": int(record["split_role"]),
                "original_cases": len(idx),
                "moved_cases": diagnostic["successful"],
                "unmoved_fallbacks": diagnostic["unmoved"],
            })
    if not np.array_equal(moved_xy[cases == 0], xy[cases == 0]):
        raise ValueError("Frozen survey background changed")
    shifts = np.concatenate(case_shifts)
    fraction = float(np.count_nonzero(shifts) / len(shifts))
    nonzero_median = float(np.median(shifts[shifts > 0])) if np.any(shifts > 0) else 0.0
    if fraction < rules["nonzero_required_fraction"] or nonzero_median < rules["global_nonzero_median_km_at_least"]:
        raise ValueError("Frozen Stage21 50-km jitter gate no longer passes")
    return moved_xy, per_species, {
        "shifted_cases": int(np.count_nonzero(shifts)),
        "fallback_cases": int(np.count_nonzero(shifts == 0)),
        "shift_median_km": float(np.median(shifts)),
        "shift_p95_km": float(np.quantile(shifts, 0.95)),
        "fraction_shifted": fraction,
    }


def jittered_predictor_support(xy, cases, sp, records, moved_xy, per_species, rules, work, climate_assets):
    """Re-evaluate exactly the original 8 climate + 5 GIS predictors, at moved sites."""
    features, names, receipts, details = [], [], [], {}
    climate = json.loads(climate_assets.read_text())
    if len(climate["assets"]) != 8:
        raise ValueError("Frozen CHELSA asset count changed")
    for i, asset in enumerate(climate["assets"]):
        path = work / f"stage22_climate_{i}.tif"
        obtain(asset["url"], path, asset["size_bytes"], asset["sha256"])
        val = raster_values(path, moved_xy)
        features.append(val)
        name = f"{asset['variable']}_{asset['time_index']}"
        names.append(name)
        receipts.append({"role": name, "bytes": asset["size_bytes"], "sha256": asset["sha256"]})
        path.unlink()
        print("JITTER_CLIMATE_VERIFIED", name, int(np.isfinite(val).sum()), flush=True)
    baseline = load("frozen_sampling_barriers_history_v01.json")
    gis = load("original_barrier_baseline_source_candidates_v01.json")
    source_sha = {a["term"]: a for a in load("barrier_gis_publisher_SHA_lock_v01.json")["asset_hashes"]}
    vectors = {}
    for src in gis["source_candidates"]:
        term = src["term"]
        receipt = source_sha[term]
        path = work / f"stage22_{term}.zip"
        obtain(src["url"], path, receipt["bytes"], receipt["sha256"])
        receipts.append({"role": term, "bytes": receipt["bytes"], "sha256": receipt["sha256"]})
        if term == "unep_mountains":
            with zipfile.ZipFile(path) as z:
                if any(n.startswith("/") or ".." in Path(n).parts for n in z.namelist()):
                    raise ValueError("Unsafe restricted GIS ZIP path")
                z.extractall(work / "mountains")
            native_grid = work / "mountains/Mountain_UNEP-WCMC_2002/Mountain_UNEP-WCMC_2002/geog/mtn_02g"
            val, detail = mountain_distances(native_grid, moved_xy)
            features.append(val)
            names.append("mountain_km")
            details[term] = detail
        else:
            vectors[term], details[term] = vector_geometries(path, term, baseline)
        path.unlink()
        print("JITTER_GIS_SOURCE_VERIFIED", term, flush=True)
    pairs = [
        ("water_km", vectors["naturalearth_rivers"] + vectors["naturalearth_lakes"]),
        ("ecoregion1_km", vectors["ecoregion_level1"]),
        ("ecoregion2_km", vectors["ecoregion_level2"]),
    ]
    for role, geoms in pairs:
        vertices = dense_vertices(geoms)
        val = vertex_distances(moved_xy, vertices)
        if role == "water_km":
            for lake in vectors["naturalearth_lakes"]:
                val[shapely.intersects_xy(lake, moved_xy[:, 0], moved_xy[:, 1])] = 0
        features.append(val)
        names.append(role)
        details[role] = {"dense_vertices": len(vertices), "max_spherical_discretization_km": 0.79}
        del vertices
        gc.collect()
        print("JITTER_BARRIER_VERIFIED", role, flush=True)
    tri = baseline["source_gis"][-1]
    tri_path = work / "stage22_tri.tif"
    obtain(tri["url"], tri_path, tri["bytes"], tri["sha256"])
    features.append(raster_values(tri_path, moved_xy))
    names.append("TRI")
    receipts.append({"role": "TRI", "bytes": tri["bytes"], "sha256": tri["sha256"]})
    tri_path.unlink()
    x = np.column_stack(features)
    if x.shape != (76700, 13):
        raise ValueError("Jittered 13-feature matrix dimensions changed")
    complete_climate = complete_rows(x[:, :8])
    complete_all = complete_rows(x)
    eligible = []
    for i, item in enumerate(per_species):
        sel = sp == i
        pos = sel & (cases == 1)
        neg = sel & (cases == 0)
        npos = int((pos & complete_all).sum())
        nneg = int((neg & complete_all).sum())
        item.update({
            "climate_complete_cases": int((pos & complete_climate).sum()),
            "all_features_complete_cases": npos,
            "all_features_complete_background": nneg,
            "eligible_after_jitter": npos >= 10 and nneg >= 20 and nneg >= npos,
        })
        if item["eligible_after_jitter"]:
            eligible.append(item)
    roles = [sum(item["role"] == role for item in eligible) for role in (0, 1)]
    passed = len(eligible) >= 50 and min(roles) >= 20
    return {
        "status": "PASS_JITTERED_INPUT_SUPPORT_NOT_EMPIRICAL_EFFECT" if passed else "NOT_EVALUABLE_JITTERED_INPUT_SUPPORT",
        "feature_names": names,
        "finite_by_feature": {name: int(np.isfinite(x[:, i]).sum()) for i, name in enumerate(names)},
        "all_eight_climates_complete": int(complete_climate.sum()),
        "all_thirteen_features_complete": int(complete_all.sum()),
        "eligible_species": len(eligible),
        "eligible_species_roles": roles,
        "source_receipts": receipts,
        "barrier_details": details,
        "per_species_aggregate": per_species,
        "moved_coordinate_sha256": sha_bytes(moved_xy, "<f8"),
        "moved_feature_matrix_sha256": sha_bytes(x, "<f8"),
    }


def run(args):
    args.work.mkdir(parents=True, exist_ok=True)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    panel = json.loads(args.panel.read_text())
    frozen = load("stage21_frozen_spatial_calibration_handoff_v01.json")
    if frozen["status"] != "FROZEN_PASS_SPATIAL_NULL_CALIBRATION_AND_GEOMETRIC_JITTER":
        raise ValueError("Stage21 source-calibration authorization missing")
    jitter_contract = load("stage21_jitter_geometry_contract_v01.json")
    xy, cases, sp, records = rebuild_sites(panel, args.work)
    if args.stage20_report is not None:
        old = json.loads(args.stage20_report.read_text())
        original = {s["species"]: s for s in old["species_support"]}
        drifts = [dict(species=s["scientific_name"], current_sha=s["coordinate_sha256"],
                       old_sha=original.get(s["scientific_name"], {}).get("coordinate_sha256"))
                  for s in records if s["coordinate_sha256"] !=
                  original.get(s["scientific_name"], {}).get("coordinate_sha256")]
        print("FROZEN_PER_SPECIES_COORDINATE_DIAGNOSTIC", json.dumps(dict(
            same_blocks=len(records)-len(drifts), changed_blocks=len(drifts),
            first_ten=drifts[:10]), sort_keys=True), flush=True)
    roles = ensure_frozen_geometry(xy, cases, sp, records, panel)
    moved, species, shift = jitter_cases_from_verified_sources(
        xy, cases, sp, records, args.work, jitter_contract
    )
    support = jittered_predictor_support(
        xy, cases, sp, records, moved, species, jitter_contract,
        args.work, args.climate_assets
    )
    report = {
        "schema": "ttf_jensen80_stage22_jittered_predictor_coverage_v01",
        "base_stage21_source_handoff_sha": hashlib.sha256(
            (Path(__file__).resolve().parent / "stage21_frozen_spatial_calibration_handoff_v01.json").read_bytes()
        ).hexdigest(),
        "base_replayed_feature_sha256": REPLAY_FEATURE_SHA256,
        "frozen_original_coordinate_sha256": ORIGINAL_XY_SHA256,
        "frozen_case_label_sha256": ORIGINAL_CASE_SHA256,
        "species": 80,
        "species_disjoint_roles": roles,
        "sites": 76700,
        "cases": 10013,
        "fixed_survey_backgrounds": 66687,
        "jitter": shift,
        **support,
        "restricted_geometries_or_coordinates_exported": False,
        "real_break_AUC_opened": False,
        "genetic_model_fitted": False,
        "prior_closed_genetic_responses_reopened": False,
        "model_opening_authorized_by_this_gate_alone": False,
        "software": {
            "python": platform.python_version(),
            "numpy": np.__version__,
            "scipy": scipy.__version__,
            "rasterio": rasterio.__version__,
            "shapely": shapely.__version__,
            "pyproj": pyproj.__version__,
        },
    }
    args.output.write_text(json.dumps(report, sort_keys=True, indent=2) + "\n")
    print(json.dumps({
        k: report[k] for k in (
            "status", "eligible_species", "eligible_species_roles",
            "all_thirteen_features_complete", "jitter"
        )
    }, sort_keys=True), flush=True)
    if report["status"] != "PASS_JITTERED_INPUT_SUPPORT_NOT_EMPIRICAL_EFFECT":
        raise SystemExit("STOP: non-evaluable jittered predictor source. Do not open AUC.")


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    for arg in ("panel", "climate-assets", "work", "output"):
        p.add_argument("--" + arg, type=Path, required=True)
    p.add_argument("--stage20-report", type=Path, required=False)
    run(p.parse_args())
