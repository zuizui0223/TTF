"""Stage23: fixed 50km case-position sensitivity of LGM-to-modern predictor magnitudes.

No genetic response values, fitted models, species-heldout AUC, or leakage into
closed outcomes. Source site coordinates and pointwise climate stay only in an
ephemeral GitHub Actions runner; public output is aggregate QA and quantiles.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import platform
from pathlib import Path

import numpy as np
from scipy.stats import spearmanr

from extract_frozen_site_features import obtain
from site_feature_core import raster_values
from stage22_jittered_predictor_support import (
    ensure_frozen_geometry, jitter_cases_from_verified_sources,
    rebuild_sites, sha_bytes,
)

HERE = Path(__file__).resolve().parent
EXPECTED_FIELDS = [
    "bio01_-190", "bio07_-190", "bio12_-190", "bio15_-190",
    "bio01_20", "bio07_20", "bio12_20", "bio15_20",
]


def case_history_magnitude(climate8, frozen_current_sd):
    x = np.asarray(climate8, dtype=np.float64)
    sd = np.asarray(frozen_current_sd, dtype=np.float64)
    if x.ndim != 2 or x.shape[1] != 8 or sd.shape != (4,):
        raise ValueError("Eight paired climate bands and four source-frozen SDs required")
    if not np.isfinite(x).all() or not np.isfinite(sd).all() or (sd <= 0).any():
        raise ValueError("Finite historical/current climate values and positive frozen SDs required")
    return np.linalg.norm((x[:, :4] - x[:, 4:]) / sd, axis=1)


def stable_rho(x, y):
    if len(x) < 3 or np.max(x) - np.min(x) < 1e-12 or np.max(y) - np.min(y) < 1e-12:
        return None
    rho = float(spearmanr(x, y).statistic)
    return rho if np.isfinite(rho) else None


def pair_summary(old, new, species, records, bands, frozen_sd):
    if old.shape != new.shape or old.shape[1] != 8 or len(species) != len(old):
        raise ValueError("Paired case climate rows cannot change counts or order")
    if not np.isfinite(old).all() or not np.isfinite(new).all():
        raise ValueError("No imputation for original or moved CHELSA climate")
    ho = case_history_magnitude(old, frozen_sd)
    hn = case_history_magnitude(new, frozen_sd)
    delta = np.abs(ho - hn)
    summaries = []
    for i, s in enumerate(records):
        mask = species == i
        if int(mask.sum()) != int(s["case_points"]):
            raise ValueError("Frozen case count changed for species " + str(i))
        a, b = ho[mask], hn[mask]
        summaries.append({
            "species": s["scientific_name"],
            "role": int(s["split_role"]),
            "case_points": len(a),
            "original_historical_magnitude_median": float(np.median(a)),
            "moved_historical_magnitude_median": float(np.median(b)),
            "paired_abs_history_shift_median": float(np.median(np.abs(b-a))),
            "paired_abs_history_shift_p95": float(np.quantile(np.abs(b-a), .95)),
            "within_species_spearman": stable_rho(a, b),
        })
    rhos = [x["within_species_spearman"] for x in summaries if x["within_species_spearman"] is not None]
    byband=[]
    for k,name in enumerate(bands):
        d=np.abs(new[:,k]-old[:,k])
        byband.append({
            "band":name,
            "paired_abs_change_median":float(np.median(d)),
            "paired_abs_change_p95":float(np.quantile(d,.95)),
            "spearman_original_moved":stable_rho(old[:,k],new[:,k]),
        })
    return {
        "source_qualified_species":len(records),
        "original_case_records":len(ho),
        "all_case_bands_finite":True,
        "original_history_magnitude_median":float(np.median(ho)),
        "original_history_magnitude_iqr":[float(v) for v in np.quantile(ho,[.25,.75])],
        "moved_history_magnitude_median":float(np.median(hn)),
        "moved_history_magnitude_iqr":[float(v) for v in np.quantile(hn,[.25,.75])],
        "paired_abs_history_shift_median":float(np.median(delta)),
        "paired_abs_history_shift_p95":float(np.quantile(delta,.95)),
        "original_moved_history_spearman":stable_rho(ho,hn),
        "species_median_paired_abs_history_shift":float(np.median([x["paired_abs_history_shift_median"] for x in summaries])),
        "species_median_within_species_spearman":float(np.median(rhos)) if rhos else None,
        "species_with_undefined_spearman":len(records)-len(rhos),
        "per_band_sensitivity":byband,
        "per_species_aggregates":summaries,
    }


def run(args):
    contract = json.loads(args.contract.read_text())
    if contract.get("status") != "LOCKED_BEFORE_PAIRED_CASE_CLIMATE_SHIFT_SCORES_NO_GENETIC_AUC":
        raise ValueError("Stage23 source-only protocol not frozen")
    parent = contract["parent"]
    metrics = contract["standardization"]
    if len(metrics["frozen_training_current_sd"]) != 4:
        raise ValueError("Stage21 frozen climatic scale missing")
    args.work.mkdir(parents=True,exist_ok=True)
    args.output.parent.mkdir(parents=True,exist_ok=True)
    panel = json.loads(args.panel.read_text())
    xy,case,sp,records = rebuild_sites(panel,args.work)
    ensure_frozen_geometry(xy,case,sp,records,panel)
    if sha_bytes(xy,"<f8") != parent["frozen_stage20_original_coordinate_sha256"]:
        raise ValueError("Original frozen coordinate SHA changed")
    jitter = json.loads((HERE/"stage21_jitter_geometry_contract_v01.json").read_text())
    moved_xy,diagnostic,shift = jitter_cases_from_verified_sources(
        xy,case,sp,records,args.work,jitter
    )
    if sha_bytes(moved_xy,"<f8") != parent["frozen_stage22_moved_coordinate_sha256"]:
        raise ValueError("Stage22 frozen moved coordinates not reproduced")
    idx = np.flatnonzero(case==1)
    if len(idx)!=10013 or shift["shifted_cases"]!=10013:
        raise ValueError("Frozen case count or successful jitter changed")
    original_case_xy = xy[idx]
    moved_case_xy = moved_xy[idx]
    climate = json.loads(args.climate_assets.read_text())
    assets = climate["assets"]
    names = [f"{a['variable']}_{a['time_index']}" for a in assets]
    if names != EXPECTED_FIELDS:
        raise ValueError("Original 8-band CHELSA source ordering changed")
    orig_fields, moved_fields, receipts = [],[],[]
    for i,a in enumerate(assets):
        p = args.work / f"stage23_source_{i}.tif"
        obtain(a["url"],p,a["size_bytes"],a["sha256"])
        orig = raster_values(p, original_case_xy)
        mov = raster_values(p, moved_case_xy)
        orig_fields.append(orig)
        moved_fields.append(mov)
        receipts.append({"field":names[i],"bytes":a["size_bytes"],"sha256":a["sha256"],
                         "original_finite":int(np.isfinite(orig).sum()),"moved_finite":int(np.isfinite(mov).sum())})
        p.unlink()
        print("STAGE23_PAIRED_CLIMATE_SOURCE_OK",names[i],
              int(np.isfinite(orig).sum()),int(np.isfinite(mov).sum()),flush=True)
    o=np.column_stack(orig_fields)
    m=np.column_stack(moved_fields)
    s=sp[idx]
    complete=bool(np.isfinite(o).all() and np.isfinite(m).all())
    output={
        "schema":"ttf_jensen80_stage23_original_vs_50km_jittered_case_climate_source_v01",
        "contract_status":contract["status"],
        "frozen_original_coordinate_sha256":sha_bytes(xy,"<f8"),
        "frozen_moved_coordinate_sha256":sha_bytes(moved_xy,"<f8"),
        "frozen_case_sha256":hashlib.sha256(case.tobytes()).hexdigest(),
        "original_case_climate_sha256":sha_bytes(o,"<f8"),
        "moved_case_climate_sha256":sha_bytes(m,"<f8"),
        "species":80,
        "training_species":44,
        "heldout_species":36,
        "case_site_records":len(idx),
        "frozen_background_records":int((case==0).sum()),
        "jitter_geometry_unchanged":shift,
        "source_receipts":receipts,
        "status":"PASS_PAIRED_CLIMATE_PREDICTOR_SENSITIVITY_MEASURED" if complete else "NOT_EVALUABLE_PAIRED_CASE_CLIMATE",
        "genetic_break_AUC_opened":False,
        "genetic_model_fitted":False,
        "pointwise_cases_or_climate_exported":False,
        "source_reconstructed_gis_not_Jensen_original_csv":True,
        "supports_historical_genetic_effect_claim":False,
        "software":{"python":platform.python_version(),"numpy":np.__version__},
    }
    if complete:
        output["summary"]=pair_summary(o,m,s,records,names,metrics["frozen_training_current_sd"])
    else:
        output["missing_pairs_total"] = int(np.count_nonzero(~np.isfinite(o).all(axis=1) | ~np.isfinite(m).all(axis=1)))
    args.output.write_text(json.dumps(output,indent=2,sort_keys=True)+"\n")
    print(json.dumps({k:output[k] for k in
                     ("status","species","case_site_records","original_case_climate_sha256","moved_case_climate_sha256")},
                     sort_keys=True),flush=True)
    if not complete:
        raise SystemExit("STOP before genetic outcome: historical-climate site sensitivity not evaluable")


if __name__=="__main__":
    p=argparse.ArgumentParser()
    for a in ("panel","climate-assets","contract","work","output"):
        p.add_argument("--"+a,type=Path,required=True)
    run(p.parse_args())
