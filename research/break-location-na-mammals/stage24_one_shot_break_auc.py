"""First and only qualified empirical Jensen80 break-LOCATION species-heldout AUC.

The real response is the original mapped break-site case vs the pre-frozen
surveyed-range background. Input SHA, original cohort, source-geography split,
synthetic known-truth qualification, and independent jitter numeric replay
MUST pass BEFORE either logistic model is fitted.

Every original/IUCN-derived site grid and pointwise source feature is ephemeral.
Never upload or log pointwise data or restricted underlying source material.
"""
from __future__ import annotations

import argparse
import gc
import hashlib
import json
import platform
from pathlib import Path
import zipfile

import numpy as np
import pyproj
import rasterio
import scipy
import shapely
from scipy.spatial import cKDTree

from extract_frozen_site_features import (
    obtain, load, vector_geometries, mountain_distances, rebuild_sites
)
from site_feature_core import dense_vertices, raster_values, vertex_distances, complete_rows
from stage22_jittered_predictor_support import (
    ensure_frozen_geometry, jitter_cases_from_verified_sources, sha_bytes
)
from spatial_null_core import (
    BASELINE_COEFF,standardize_predictors,training_history_unique_fraction,
    evaluate_species_holdout
)

HERE=Path(__file__).resolve().parent
NAMES=["bio01_-190","bio07_-190","bio12_-190","bio15_-190",
       "bio01_20","bio07_20","bio12_20","bio15_20",
       "mountain_km","water_km","ecoregion1_km","ecoregion2_km","TRI"]


def mandatory_authorities(contract):
    if contract.get("status")!="FROZEN_PRE_EMPIRICAL_GENETIC_BREAK_LOCATION_AUC":
        raise ValueError("Source/model original one-shot not frozen")
    auth=contract["authority"]
    s21=load("stage21_frozen_spatial_calibration_handoff_v01.json")
    if s21.get("status")!="FROZEN_PASS_SPATIAL_NULL_CALIBRATION_AND_GEOMETRIC_JITTER":
        raise ValueError("Synthetic spatial known-truth gate not passed")
    synthetic=s21["calibration"]
    if (synthetic["null_false_positives"],synthetic["n_null"],
        synthetic["positive_detections"],synthetic["n_positive"])!=(5,200,200,200):
        raise ValueError("Frozen synthetic family mismatch")
    s22=load("stage22_two_run_numeric_replay_amendment_v01.json")
    if s22.get("status")!="TWO_INDEPENDENT_REPLAYS_BYTE_EQUAL_FIRST_SUCCESS_NUMERIC_SHA_DISCREPANCY_PRESERVED":
        raise ValueError("Source-only two-run jitter numeric replay not yet accepted")
    for field in ("exact_replay_a","exact_replay_b"):
        if s22[field]["moved_feature_matrix_sha256"]!=auth["jittered_expected_numeric_matrix_sha256"]:
            raise ValueError("Independent jittered numerical replay SHA not bound")
    s23=load("stage23_case_climate_sensitivity_result_v01.json")
    if s23.get("status")!="PASS_PAIRED_CLIMATE_PREDICTOR_SENSITIVITY_MEASURED_NOT_GENETIC_EFFECT":
        raise ValueError("Original vs 50km climate predictor sensitivity not audited")
    s21c=load("stage21_spatial_null_contract_v01.json")
    if s21c["status"]!="FROZEN_BEFORE_ANY_GENETIC_HISTORY_AUC":
        raise ValueError("Frozen Stage21 model settings changed")
    if not np.array_equal(BASELINE_COEFF,np.asarray(s21c["baseline_coefficients_in_training_SD_units"])):
        raise ValueError("Synthetic model source differs")
    if contract["model"]["bootstrap"].split(";")[0]!="1,999 species cluster replicates":
        raise ValueError("Original bootstrap rule changed")
    return s21c


def original_and_jittered_features(xy,moved_xy,work,assets):
    """Same math/source code as frozen Stage20/22; separate queries preserve exact raw SHAs."""
    orig,moved,names,receipts=[],[],[],[]
    climate=json.loads(Path(assets).read_text())
    if len(climate["assets"])!=8:
        raise ValueError("Eight and only eight frozen original CHELSA fields required")
    for i,a in enumerate(climate["assets"]):
        p=work/f"stage24_climate_{i}.tif"
        obtain(a["url"],p,a["size_bytes"],a["sha256"])
        orig.append(raster_values(p,xy))
        moved.append(raster_values(p,moved_xy))
        name=f"{a['variable']}_{a['time_index']}"
        names.append(name)
        receipts.append({"role":name,"size_bytes":a["size_bytes"],"sha256":a["sha256"]})
        p.unlink()
        print("PRE_RESPONSE_CHELSEA_SOURCE_VERIFIED",name,flush=True)
    contract=load("frozen_sampling_barriers_history_v01.json")
    gis=load("original_barrier_baseline_source_candidates_v01.json")
    sha={a["term"]:a for a in load("barrier_gis_publisher_SHA_lock_v01.json")["asset_hashes"]}
    vectors={}
    for src in gis["source_candidates"]:
        term=src["term"];h=sha[term];p=work/(term+".zip")
        obtain(src["url"],p,h["bytes"],h["sha256"])
        receipts.append({"role":term,"size_bytes":h["bytes"],"sha256":h["sha256"]})
        if term=="unep_mountains":
            with zipfile.ZipFile(p) as z:
                if any(n.startswith("/") or ".." in Path(n).parts for n in z.namelist()):
                    raise ValueError("Unsafe source ZIP path")
                z.extractall(work/"mountains")
            grid=work/"mountains/Mountain_UNEP-WCMC_2002/Mountain_UNEP-WCMC_2002/geog/mtn_02g"
            a,_=mountain_distances(grid,xy)
            b,_=mountain_distances(grid,moved_xy)
            orig.append(a);moved.append(b);names.append("mountain_km")
        else:
            vectors[term],_=vector_geometries(p,term,contract)
        p.unlink()
        print("PRE_RESPONSE_GIS_SOURCE_VERIFIED",term,flush=True)
    for name,geom in [
        ("water_km",vectors["naturalearth_rivers"]+vectors["naturalearth_lakes"]),
        ("ecoregion1_km",vectors["ecoregion_level1"]),
        ("ecoregion2_km",vectors["ecoregion_level2"])
    ]:
        vertices=dense_vertices(geom)
        a=vertex_distances(xy,vertices)
        b=vertex_distances(moved_xy,vertices)
        if name=="water_km":
            for lake in vectors["naturalearth_lakes"]:
                a[shapely.intersects_xy(lake,xy[:,0],xy[:,1])]=0
                b[shapely.intersects_xy(lake,moved_xy[:,0],moved_xy[:,1])]=0
        orig.append(a);moved.append(b);names.append(name)
        del vertices
        gc.collect()
        print("PRE_RESPONSE_BARRIER_VALUES_VERIFIED",name,flush=True)
    tri=contract["source_gis"][-1]
    p=work/"stage24_tri.tif"
    obtain(tri["url"],p,tri["bytes"],tri["sha256"])
    orig.append(raster_values(p,xy))
    moved.append(raster_values(p,moved_xy))
    names.append("TRI")
    receipts.append({"role":"TRI","size_bytes":tri["bytes"],"sha256":tri["sha256"]})
    p.unlink()
    x0=np.column_stack(orig)
    x1=np.column_stack(moved)
    if names!=NAMES or x0.shape!=(76700,13) or x1.shape!=(76700,13):
        raise ValueError("The frozen 13-column source order or source row geometry changed")
    if not complete_rows(x0).all() or not complete_rows(x1).all():
        raise ValueError("Source feature fields non-finite; no post-result imputation allowed")
    if np.any(x0[:,8:]<0) or np.any(x1[:,8:]<0):
        raise ValueError("Geographic barrier/ruggedness values invalid")
    return x0,x1,receipts


def frozen_original_predictor_scale_for_jitter(xm,lonlat,mu,scale,climate_sd):
    """Apply EXACT original training species standardization to jittered features."""
    h=np.linalg.norm((xm[:,:4]-xm[:,4:8])/climate_sd,axis=1)
    p=np.column_stack([xm[:,4:8],np.sqrt(xm[:,8:12]),np.log1p(xm[:,12]),
                       lonlat[:,1],lonlat[:,0],h])
    return (p-mu)/scale


def run(args):
    contract=json.loads(args.contract.read_text())
    sc=mandatory_authorities(contract)
    args.work.mkdir(parents=True,exist_ok=True)
    args.output.parent.mkdir(parents=True,exist_ok=True)
    panel=json.loads(args.panel.read_text())
    xy,y,sp,records=rebuild_sites(panel,args.work)
    roles=ensure_frozen_geometry(xy,y,sp,records,panel)
    if roles!=[44,36]:
        raise ValueError("Fixed species source roles changed")
    jitter=load("stage21_jitter_geometry_contract_v01.json")
    moved_xy,_,moved_count=jitter_cases_from_verified_sources(
        xy,y,sp,records,args.work,jitter
    )
    auth=contract["authority"]
    if sha_bytes(xy,"<f8")!=auth["original_expected_coord_sha256"]:
        raise ValueError("Original frozen coordinate SHA not reproduced")
    if sha_bytes(moved_xy,"<f8")!=auth["jittered_expected_coord_sha256"]:
        raise ValueError("Jittered frozen coordinate SHA not reproduced")
    if hashlib.sha256(y.tobytes()).hexdigest()!=auth["case_labels_expected_sha256"]:
        raise ValueError("Frozen case labels do not reproduce")
    if len(xy)!=contract["source_geometry"]["total"] or int(y.sum())!=10013:
        raise ValueError("Case/background source counts changed")
    x0,x1,receipts=original_and_jittered_features(
        xy,moved_xy,args.work,args.climate_assets
    )
    src0=sha_bytes(x0,"<f8")
    src1=sha_bytes(x1,"<f8")
    # Critical FIREWALL: Do not look at any genetic-response model unless
    # BOTH original+50km source numeric matrices match prospectively pinned
    # independent SHA-verified source-replay conditions.
    if src0!=auth["original_expected_numeric_matrix_sha256"]:
        raise ValueError("STOP BEFORE BIOLOGICAL RESPONSE: original 13-field source numeric SHA differs: "+src0)
    if src1!=auth["jittered_expected_numeric_matrix_sha256"]:
        raise ValueError("STOP BEFORE BIOLOGICAL RESPONSE: jittered 13-field source numeric SHA differs: "+src1)
    if len(receipts)!=14:
        raise ValueError("Not all 14 original climate/geography provenance sources verified")
    role_by_species=np.asarray([int(s["split_role"]) for s in records],dtype=int)
    role=role_by_species[sp]
    if not np.array_equal(np.unique(sp),np.arange(80)):
        raise ValueError("Source species identifiers drifted")
    p0,mu,scale,climate_sd=standardize_predictors(x0,xy,sp,role)
    uniq=training_history_unique_fraction(p0,sp,role)
    expected_unique=load("stage21_frozen_spatial_calibration_handoff_v01.json")["calibration"]["training_history_unique_fraction"]
    if not abs(uniq-expected_unique)<1e-8:
        raise ValueError("STOP BEFORE BIOLOGICAL RESPONSE: synthetic unique predictor fraction changed")
    p1=frozen_original_predictor_scale_for_jitter(x1,moved_xy,mu,scale,climate_sd)
    if not np.isfinite(p1).all():
        raise ValueError("Jitter predictor nonfinite after frozen scaling")
    # Every source/data qualification passed. The following is the ONE and
    # ONLY authorized response opening for the entirely distinct mapped
    # genetic break-location endpoint. Never rerun after an empirical result.
    print("AUTHORIZED_ONE_SHOT_REAL_BREAK_LOCATION_MODEL_BEGIN",flush=True)
    primary=evaluate_species_holdout(p0,y,sp,role,bootstrap_reps=1999,seed=20261029)
    sensitivity=evaluate_species_holdout(p1,y,sp,role,bootstrap_reps=1999,seed=20261029)
    a=primary["positive_decision"]
    b=sensitivity["positive_decision"]
    decision=("QUALIFIED_POSITIVE_HISTORY_BREAK_LOCATION_PREDICTION_ASSOCIATION_WITH_50KM_SENSITIVITY"
              if a and b else
              "POSITION_SENSITIVE_NOT_ROBUST_HISTORY_BREAK_LOCATION_ASSOCIATION" if a
              else "NO_QUALIFIED_ADDED_HISTORY_LOCATION_PREDICTION")
    report={
        "schema":"ttf_jensen80_stage24_first_authorized_empirical_genetic_break_location_AUC_v01",
        "status":decision,
        "empirical_outcome_opened_once":True,
        "source_model_contract":"stage24_empirical_break_location_one_shot_contract_v01.json",
        "n_species":80,"train_species":44,"heldout_species":36,
        "case_site_records":int(y.sum()),"background_records":int((y==0).sum()),
        "input_original_feature_sha256":src0,
        "input_jittered_feature_sha256":src1,
        "input_original_coordinate_sha256":sha_bytes(xy,"<f8"),
        "input_jittered_coordinate_sha256":sha_bytes(moved_xy,"<f8"),
        "source_receipts":receipts,
        "training_history_unique_fraction":uniq,
        "prior_synthetic_null_positive_worlds_unchanged":True,
        "original_case_location_primary":primary,
        "frozen_50km_case_location_sensitivity":sensitivity,
        "input_positional_jitter_km":moved_count,
        "per_program_positive_requires_9833_species_cluster_CI":True,
        "fitted_source_Jensen_original_author_baseline_exactly":False,
        "backgrounds_are_availability_not_true_break_absence":True,
        "causal_LGM_barrier_refugia_effect_identified":False,
        "prior_He_Fst_COI_TTFQ_endpoints_reopened":False,
        "licensed_IUCN_UNEP_or_pointwise_output_exported":False,
        "software":{"python":platform.python_version(),"numpy":np.__version__,
                    "scipy":scipy.__version__,"rasterio":rasterio.__version__,
                    "shapely":shapely.__version__,"pyproj":pyproj.__version__}
    }
    args.output.write_text(json.dumps(report,sort_keys=True,indent=2)+"\n")
    print(json.dumps({k:report[k] for k in ("status","original_case_location_primary",
                                        "frozen_50km_case_location_sensitivity")},sort_keys=True),flush=True)


if __name__=="__main__":
    p=argparse.ArgumentParser()
    for arg in ("panel","climate-assets","contract","work","output"):
        p.add_argument("--"+arg,type=Path,required=True)
    run(p.parse_args())
