"""Response-blind information added by 21ka-to-present climate in 70 plants.

Strict input identity: same 70 source-qualified species, 104 studies, original
source's 1561 geolocated comparison records, 253 unique coordinates, frozen
35/35 species-and-study split, exact original 8 CHELSA numeric matrix SHA.

Use unique study–species–ordered-site-pair predictor geometry, equal weight
per training species. No field of plant trait change, measured rate, mean,
variance, genotype or fitness is accessed. No biological response opened.
"""
from __future__ import annotations

import argparse
from collections import Counter, defaultdict
import hashlib
import json
from pathlib import Path

import numpy as np

from proceed_v62_metadata_row_census import fetch_original_datafile_bytes
from proceed_v62_borealis_header_only import validated_get_catalog, verify_frozen_catalog
from proceed_v62_genetic_geometry_intersection import SOURCE_CONTRACT, NATIVE_LOCK
from proceed_plant_chelsa_source_support import (
    source_geometry, component_assignment, verified_raster, source_climate_coverage,
)
from site_feature_core import raster_values

HERE = Path(__file__).resolve().parent
CONTRACT = HERE / "proceed70_plant_history_information_contract_v01.json"
ASSETS = HERE / "chelsa_trace21k_fixed8_assets_v01.json"
EARTH_RADIUS_KM = 6371.0088


def deduplicate_study_species_site_pairs(cohort):
    """One source predictor geometry per exact pair, not one per trait/effect."""
    return sorted(set((sp, study, a, b) for sp, study, a, b in cohort))


def haversine_pairs(a, b):
    """Geodesic great-circle proximity; no raw pairwise coordinates are output."""
    la1 = np.deg2rad(a[:, 1])
    la2 = np.deg2rad(b[:, 1])
    dlo = np.deg2rad((b[:, 0] - a[:, 0] + 180) % 360 - 180)
    dla = la2 - la1
    h = np.sin(dla / 2)**2 + np.cos(la1) * np.cos(la2) * np.sin(dlo/2)**2
    return EARTH_RADIUS_KM * 2 * np.arcsin(np.sqrt(np.clip(h, 0, 1)))


def create_pair_input(unique, coords, original8, role_by_species):
    index = {tuple(pt):i for i, pt in enumerate(coords)}
    i=np.array([index[a] for _, _, a, b in unique], int)
    j=np.array([index[b] for _, _, a, b in unique], int)
    a=original8[i]
    b=original8[j]
    z1=np.asarray([pair[2] for pair in unique],float)
    z2=np.asarray([pair[3] for pair in unique],float)
    species=[pair[0] for pair in unique]
    studies=[pair[1] for pair in unique]
    role=np.array([role_by_species[s] for s in species])
    return a,b,z1,z2,species,studies,role


def equal_species_weights(species, role):
    n=Counter(s for s,r in zip(species,role) if r=="train")
    k=len(n)
    if not k:raise ValueError("Training species missing")
    w=np.zeros(len(species),float)
    for i,(sp,r) in enumerate(zip(species,role)):
        if r=="train":
            w[i]=1.0/(k*n[sp])
    if abs(float(np.sum(w))-1)>1e-11:
        raise ValueError("Equal-species training weights not normalized")
    return w


def weighted_center_and_std(x,w):
    """Weights normalized separately by caller; columns with zero SD return 1."""
    mu=np.sum(x*w[:,None],axis=0)
    var=np.sum(((x-mu)**2)*w[:,None],axis=0)
    sd=np.sqrt(np.maximum(0,var))
    sd=np.where(sd>1e-12,sd,1.0)
    return mu,sd


def features_for_original_pairs(a,b,z1,z2,sd):
    if a.ndim!=2 or b.shape!=a.shape or a.shape[1]!=8:
        raise ValueError("Exactly 8 CHELSA source bands at both locations needed")
    if np.any(sd<=0) or np.any(~np.isfinite(sd)):
        raise ValueError("Frozen 0BP climate scales invalid")
    def one_history(x):
        return np.linalg.norm((x[:,:4]-x[:,4:8])/sd,axis=1)
    H=.5*(one_history(a)+one_history(b))
    curr1=a[:,4:8];curr2=b[:,4:8]
    midlat=.5*(z1[:,1]+z2[:,1])
    lon1=np.deg2rad(z1[:,0]);lon2=np.deg2rad(z2[:,0])
    midlon=np.rad2deg(np.arctan2(np.sin(lon1)+np.sin(lon2),
                                np.cos(lon1)+np.cos(lon2)))
    # Prespecified 11 modern and geographic controls, not source outcomes.
    controls=np.column_stack([
        (curr1+curr2)/2,
        np.abs(curr1-curr2),
        midlat,midlon,haversine_pairs(z1,z2)
    ])
    if controls.shape[1]!=11 or not np.isfinite(controls).all() or not np.isfinite(H).all():
        raise ValueError("Nonfinite original historical predictor or baseline control")
    return H,controls


def source_only_information(a,b,z1,z2,species,studies,role,minimum=0.05):
    idx=role=="train"
    if sum(idx)<36:
        raise ValueError("Need at least 36 unique source-pair train rows for 12 coefficients")
    w=equal_species_weights(species,role)
    # The source-original current climate SD is computed from BOTH original
    # endpoints with study–species location pairs rather than repeated effects.
    m1=a[:,4:8];m2=b[:,4:8]
    endpoint=np.vstack([m1,m2])
    endpoint_w=np.concatenate([w/2,w/2])
    _,sd=weighted_center_and_std(endpoint,endpoint_w)
    H,X=features_for_original_pairs(a,b,z1,z2,sd)
    if np.linalg.matrix_rank(np.column_stack([np.ones(idx.sum()),X[idx]]))<12:
        raise ValueError("Modern geography baseline rank deficient on source train pairs")
    mean,std=weighted_center_and_std(X[idx],w[idx])
    Z=(X-mean)/std
    design=np.column_stack([np.ones(len(X)),Z])
    sw=np.sqrt(w[idx])
    beta=np.linalg.lstsq(design[idx]*sw[:,None],H[idx]*sw,rcond=1e-12)[0]
    fitted=design@beta
    mu_H=float(np.sum(w[idx]*H[idx]))
    denom=float(np.sum(w[idx]*(H[idx]-mu_H)**2))
    if denom<=1e-12:
        raise ValueError("Historical displacement magnitude has no training variance")
    residual=H-fitted
    numerator=float(np.sum(w[idx]*residual[idx]**2))
    unique=float(np.clip(numerator/denom,0,1))
    residual_by_species=defaultdict(float)
    for sp,is_train,wi,ri in zip(species,idx,w,residual):
        if is_train:
            residual_by_species[sp]+=wi*float(ri*ri)
    shares=np.array(list(residual_by_species.values()),float)
    if shares.sum()<=1e-12:
        effective=0.0;largest=0.0;active=0
    else:
        shares=shares/shares.sum()
        effective=float(1/np.sum(shares**2))
        largest=float(shares.max())
        active=int((shares>1e-8).sum())
    # Heldout baseline is deliberately NOT fitted anew; direct transfer of
    # the training regression, no biological response, no model optimization.
    hold=role=="heldout"
    nh=Counter(sp for sp,ok in zip(species,hold) if ok)
    wh=np.array([1/(len(nh)*nh[s]) if ishold else 0
                 for s,ishold in zip(species,hold)])
    hm=float(np.sum(wh[hold]*H[hold]))
    hd=float(np.sum(wh[hold]*(H[hold]-hm)**2))
    holdout_proj_error_ratio=(float(np.sum(wh[hold]*residual[hold]**2)/hd)
                                if hd>1e-12 else None)
    return {
        "train_unique_historical_predictor_variance_fraction":unique,
        "trained_history_var_weighted":denom,
        "train_species_effective_residual_signal":effective,
        "train_species_max_signal_share":largest,
        "train_species_with_nonzero_residual_signal":active,
        "train_baseline_design_rank":int(np.linalg.matrix_rank(design[idx])),
        "train_weighted_current_0bp_SD":sd.tolist(),
        "heldout_training_projection_error_ratio_not_unique_fraction":holdout_proj_error_ratio,
        "train_pair_count":int(idx.sum()),
        "heldout_pair_count":int(hold.sum()),
        "train_species_count":len(set(sp for sp,ok in zip(species,idx) if ok)),
        "heldout_species_count":len(nh),
        "train_study_count":len(set(s for s,ok in zip(studies,idx) if ok)),
        "heldout_study_count":len(set(s for s,ok in zip(studies,hold) if ok)),
        "history_predictor_sample_nonzero_train_species_count":len(set(
            sp for sp,v,ok in zip(species,H,idx) if ok and abs(v)>1e-12)),
        "history_predictor_sample_nonzero_heldout_species_count":len(set(
            sp for sp,v,ok in zip(species,H,hold) if ok and abs(v)>1e-12)),
        "train_condition_met_predeclared":bool(
            unique>=minimum and active>=20 and
            len(set(sp for sp,v,ok in zip(species,H,hold) if ok and abs(v)>1e-12))>=20),
        "no_biological_response_used":True,
    }


def run(args):
    contract=json.loads(CONTRACT.read_text())
    if contract["status"]!="FROZEN_BEFORE_PLANT_HISTORY_INFORMATION_RESPONSE_BLIND_CALCULATION":
        raise ValueError("Prospective history information protocol unfrozen")
    locked=contract["authoritative_input"]
    src=json.loads(SOURCE_CONTRACT.read_text())["current_source"]
    original=json.loads(NATIVE_LOCK.read_text())["native_original_transport"]
    args.work.mkdir(parents=True,exist_ok=True)
    args.output.parent.mkdir(parents=True,exist_ok=True)
    report={
        "schema":"ttf_proceed70_response_blind_history_information_v01",
        "status":"NOT_EVALUABLE_SOURCE_REPLAY_CHECK_NOT_PASSED",
        "rate_or_trait_or_genetic_response_values_opened":False,
        "TTF_closed_outcomes_reopened":False,
        "raw_spatial_climate_or_response_rows_exported":False,
    }
    try:
        catalog=validated_get_catalog(
            "https://borealisdata.ca/api/datasets/:persistentId/?persistentId=doi%3A10.5683%2FSP3%2FNXSL3Q")
        verify_frozen_catalog(catalog,src)
        raw,md5,sha=fetch_original_datafile_bytes(1024206)
        if (md5!=locked["original_proceed_md5"] or sha!=locked["original_proceed_sha256"] or
            len(raw)!=original["observed_size_bytes"]):
            raise ValueError("Frozen original PROCEED source identity not reproducible")
        cohort,study_species=source_geometry(raw,contract)
        roles,split,sp_sha=component_assignment(study_species)
        if (len(cohort)!=locked["original_original_endpoint_pair_rows"] or
            len(study_species)!=locked["original_source_studies"] or
            len(roles)!=locked["original_source_plant_taxa"] or
            sp_sha!=locked["source_study_species_role_sha256"]):
            raise ValueError("Original study-species-geometry split changed")
        original_unique=deduplicate_study_species_site_pairs(cohort)
        coords=sorted({xy for _,_,a,b in cohort for xy in (a,b)})
        xy=np.asarray(coords,dtype=float)
        xy_sha=hashlib.sha256(xy.astype("<f8").tobytes()).hexdigest()
        if xy_sha!=locked["origin_lonlat_sha256"]:
            raise ValueError("Source plant site coordinate order hash changed")
        asset=json.loads(ASSETS.read_text())
        names=[f"{a['variable']}_{a['time_index']}" for a in asset["assets"]]
        if names!=["bio01_-190","bio07_-190","bio12_-190","bio15_-190",
                  "bio01_20","bio07_20","bio12_20","bio15_20"]:
            raise ValueError("Original CHELSA historical source changed")
        climate=[];receipts=[]
        for k,a in enumerate(asset["assets"]):
            path=args.work/f"p70_history_{k}.tif"
            verified_raster(a,path)
            val=raster_values(path,xy)
            climate.append(val)
            receipts.append({
                "name":names[k],"sha256":a["sha256"],"source_bytes":a["size_bytes"],
                "finite_original_sites":int(np.isfinite(val).sum()),
            })
            path.unlink()
            print("PRE_RESPONSE_PLANT_HISTORY_SOURCE_OK",names[k],flush=True)
        x=np.column_stack(climate)
        xsha=hashlib.sha256(x.astype("<f8").tobytes()).hexdigest()
        if xsha!=locked["chelsa8_complete_raw_sha256"] or not np.isfinite(x).all():
            raise ValueError("Original 8-CHELSA plant numeric matrix SHA changed or source nodata")
        a,b,z1,z2,species,studies,role=create_pair_input(original_unique,coords,x,roles)
        result=source_only_information(a,b,z1,z2,species,studies,role,minimum=
                                       contract["pass_criterion"]["minimum_training_total_unique_fraction"])
        passed=(result["train_condition_met_predeclared"] and
                result["train_pair_count"]>=36 and result["heldout_pair_count"]>0)
        report.update(result)
        report.update({
            "status":("PASS_PLANT70_HISTORY_NONREDUNDANT_INFORMATION_SOURCE_ONLY"
                      if passed else "NOT_EVALUABLE_PLANT70_HISTORICAL_RELATION_SOURCE_ONLY"),
            "original_effect_row_count":len(cohort),
            "unique_study_species_site_pair_predictor_rows":len(original_unique),
            "original_source_unique_site_count":len(coords),
            "source_study_species_split_sha256":sp_sha,
            "original_lonlat_sha256":xy_sha,
            "original_chelsa8_matrix_sha256":xsha,
            "asset_sha_receipts":receipts,
        })
    except Exception as exc:
        report["status"]="NOT_EVALUABLE_SOURCE_REPLAY_OR_HISTORY_GEOMETRY"
        report["error_type"]=type(exc).__name__
        report["error_message"]=str(exc)[:400]
    args.output.write_text(json.dumps(report,sort_keys=True,indent=2)+"\n")
    print(json.dumps({k:report.get(k) for k in (
        "status","original_effect_row_count",
        "unique_study_species_site_pair_predictor_rows",
        "train_unique_historical_predictor_variance_fraction",
        "train_species_effective_residual_signal","train_pair_count","heldout_pair_count")},
        sort_keys=True))


if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--work",required=True,type=Path)
    p.add_argument("--output",required=True,type=Path)
    run(p.parse_args())
