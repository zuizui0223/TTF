"""Prospective Haldane source-only CHELSA geometry and synthetic null/power gate.

This stage never reads, parses, calculates, or fits any original phenotype
mean/SD/sample-size magnitude, Haldane/Darwin rate, genotype or adaptation
outcome. Binary field-presence masks only, as previously frozen.
"""
from __future__ import annotations
import argparse,csv,hashlib,io,json,math
from collections import Counter,defaultdict
from pathlib import Path

import numpy as np
from proceed_v62_borealis_header_only import validated_get_catalog,verify_frozen_catalog
from proceed_v62_metadata_row_census import (fetch_original_datafile_bytes,
                                              valid_binomial,valid_coordinates)
from proceed70_haldane_presence_mask import FIELDS,META,missingness_summary,is_present_raw_token,CONTRACT as PRESENCE
from proceed_plant_chelsa_source_support import (source_geometry,component_assignment,
     verified_raster)
from proceed70_plant_history_information import (
    create_pair_input,features_for_original_pairs,source_only_information,
    weighted_center_and_std
)
from proceed70_haldane_rate_model_core import (
    equal_species_weights,scale_train,ridge_predict,species_mse_comparison,
    known_truth_labels,wilson95,SEED,ALPHA
)
from proceed_v62_genetic_geometry_intersection import SOURCE_CONTRACT,NATIVE_LOCK
from site_feature_core import raster_values

HERE=Path(__file__).resolve().parent
CONTRACT=HERE/"proceed70_haldane_first_shot_contract_v01.json"
CLIMATE=HERE/"chelsa_trace21k_fixed8_assets_v01.json"
FULL_GEO=HERE/"proceed_plant_chelsa_source_support_contract_v01.json"


def presence_keyed_metadata(raw,presence_contract):
    empty=set(presence_contract["field_presence_only"]["absent_if_trimmed_lowercase_in"])
    r=csv.reader(io.TextIOWrapper(io.BytesIO(raw),encoding="utf-8-sig",newline=""),
                 delimiter=",",strict=True)
    cols=next(r)
    meta=list(META)+["generations","experimental"]
    must=list(FIELDS)
    if not set(meta+must).issubset(cols) or len(cols)!=69:
        raise ValueError("Original PROCEED native metadata field order drifted")
    ix={name:cols.index(name) for name in set(meta+must)}
    info=defaultdict(lambda:{"generations":[],"trait":Counter(),
                             "design":Counter(),"experimental":Counter(),
                             "record_count":0})
    for raw_row in r:
        if not raw_row:continue
        if len(raw_row)!=69:raise ValueError("Native response source row width changed")
        # Metadata-only source, response fields inspected for BINARY presence:
        # Do not convert or compare a single mean/SD/sample-size value.
        m={k:raw_row[ix[k]].strip() for k in meta}
        if m["phylum"]!="Streptophyta" or m["taxa"]!="Plant" or m["genphen"]!="Genetic":
            continue
        a=valid_coordinates(m["sample1_latitude"],m["sample1_longitude"])
        b=valid_coordinates(m["sample2_latitude"],m["sample2_longitude"])
        species=valid_binomial(m["sp_pub"]) or valid_binomial(m["sp_ncbi"])
        if not(a and b and species and m["ref_ID"]):continue
        flags=[is_present_raw_token(raw_row[ix[name]],empty) for name in FIELDS]
        if not all(flags):continue
        # Generations is *environmental/time metadata*, not the trait outcome.
        try:g=float(m["generations"])
        except ValueError:raise ValueError("Non-numeric generation-time source metadata before biological response")
        if not math.isfinite(g) or g<=0:raise ValueError("Nonpositive generation count in H-presence-only source geometry")
        key=(species,m["ref_ID"],(a[1],a[0]),(b[1],b[0]))
        slot=info[key]
        slot["generations"].append(g)
        slot["trait"][m["trait_type"] or "<missing>"]+=1
        slot["design"][m["design"] or "<missing>"]+=1
        slot["experimental"][m["experimental"] or "<missing>"]+=1
        slot["record_count"]+=1
    return info


def method_features(keys,meta):
    families=("growth","otherLH","othermorphology","phenology","physio","response","size")
    out=[]
    for key in keys:
        obj=meta[key];n=obj["record_count"]
        if not n or set(obj["trait"])-set(families):
            raise ValueError("Unexpected source trait category or no admissible rows")
        if set(obj["design"])-{"Allochronic","Synchronic"} or set(obj["experimental"])-{"Y","N"}:
            raise ValueError("Unexpected study design or experimental source code")
        out.append([
          obj["design"]["Synchronic"]/n,
          obj["experimental"]["Y"]/n,
          np.log1p(np.median(obj["generations"])),
          np.log1p(n),
          *[obj["trait"][family]/n for family in families[:6]]
        ])
    a=np.asarray(out,float)
    if a.shape!=(len(keys),10) or not np.isfinite(a).all():
        raise ValueError("Frozen ten independent original-study metadata controls missing")
    return a


def source_history_concentration(H,X,sp,role):
    weights=equal_species_weights(sp,role)
    t=role=="train"
    w=weights[t]
    if t.sum()<36:raise ValueError("Insufficient H-mask source pair support")
    mu,sd=scale_train(X[t],w)
    Z=(X-mu)/sd
    design=np.column_stack([np.ones(len(H)),Z])
    if not np.isfinite(H).all():raise ValueError("Historical source predictor non-finite")
    beta=np.linalg.lstsq(design[t]*np.sqrt(w)[:,None],
                         H[t]*np.sqrt(w),rcond=1e-12)[0]
    resid=H-design@beta
    mh=w@H[t]
    total=w@((H[t]-mh)**2)
    unique=(w@(resid[t]**2))/total if total>1e-12 else 0.
    support=defaultdict(float)
    for species,ok,weight,error in zip(sp,t,weights,resid):
        if ok:support[species]+=weight*error*error
    shares=np.asarray(list(support.values()))
    if shares.sum()<=1e-12:
        eff=0;peak=0
    else:
        shares=shares/shares.sum()
        eff=float(1/np.sum(shares**2))
        peak=float(shares.max())
    return {"history_unique_fraction_after_all_21_controls":float(unique),
            "effective_source_species_of_history_residual":eff,
            "maximum_one_source_species_history_residual_share":peak,
            "source_baseline_rank_including_intercept":int(np.linalg.matrix_rank(design[t])),
            "source_baseline_dim_including_intercept":int(design.shape[1])}


def source_only_geometry(raw,climate_matrix,coords,source_contract):
    full_cohort,study_species=source_geometry(raw,source_contract)
    roles,split,digest=component_assignment(study_species)
    if len(full_cohort)!=1561 or len(roles)!=70 or digest!="73b16a46f4458aee64d735a495273e5d760683b020c5432a9701fc14214f231d":
        raise ValueError("Original pre-response 70-plant source graph changed")
    p=json.loads(PRESENCE.read_text())
    old=missingness_summary(raw,p)
    if (old["status"]!="HALDANE_PRESENCE_MASK_QUALIFIED_BUT_RATE_VALUES_NOT_OPENED"
        or old["distinct_species_with_any_complete_mask"]!=56
        or old["rows_with_all_required_Haldane_source_fields_present"]!=1279
        or old["distinct_original_study_species_geographic_pairs_with_any_complete_mask"]!=165
        or old["complete_mask_original_species_by_role"]!={"train":29,"heldout":27}):
        raise ValueError("Original 56-plant Haldane binary presence panel changed")
    meta=presence_keyed_metadata(raw,p)
    keys=sorted(meta)
    if len(keys)!=165:
        raise ValueError("Source-only Haldane geometry must have 165 unique study-species pairs")
    a,b,z1,z2,species,studies,role=create_pair_input(keys,coords,climate_matrix,roles)
    s=np.asarray(species);studies=np.asarray(studies);role=np.asarray(role)
    w=equal_species_weights(s,role)
    t=role=="train"
    # Each original site endpoint contributes half the training species weight
    endpoints=np.vstack([a[:,4:8],b[:,4:8]])
    weights=np.concatenate([w/2,w/2])
    _,sd=weighted_center_and_std(endpoints,weights)
    H,base11=features_for_original_pairs(a,b,z1,z2,sd)
    meta10=method_features(keys,meta)
    X=np.column_stack([base11,meta10])
    if X.shape!=(len(keys),21) or not np.isfinite(X).all():
        raise ValueError("Source-only 11 climate/geography + 10 methodology controls not finite")
    if len(set(s[role=="train"]))!=29 or len(set(s[role=="heldout"]))!=27:
        raise ValueError("Original 29/27 masked source roles changed")
    return X,H,s,studies,role,sd,old


def known_truth_gate(X,H,sp,studies,role,contract):
    frozen=contract["synthetic_pre_response_qualification"]
    truths={}
    for kind,beta in (("null",0.0),("positive",frozen["positive_source_history_coef_per_training_SD"])):
        yes=0;results=[]
        for i in range(frozen["simulated_worlds_"+kind]):
            seed=frozen["frozen_seed"]+i+(100000 if kind=="positive" else 0)
            rng=np.random.default_rng(seed)
            y=known_truth_labels(X,H,sp,studies,role,rng,beta)
            base=ridge_predict(X,y,sp,role)
            hist=ridge_predict(X,y,sp,role,history=H)
            decision=species_mse_comparison(y,base,hist,sp,role,
                       reps=frozen["bootstrap_draws_per_synthetic_world"],seed=seed+500000)
            yes+=bool(decision["positive"])
            results.append(decision["relative_MSE_gain"])
        n=len(results)
        truths[kind]={"worlds":n,"positive_decisions":yes,"Wilson95":wilson95(yes,n),
                     "median_heldout_MSE_gain":float(np.median(results))}
        print("KNOWN_TRUTH_SYNTHETIC_ONLY",kind,yes,n,flush=True)
    passed=(truths["null"]["Wilson95"][1]<=frozen["null_gate_Wilson95_upper_max"]
             and truths["positive"]["Wilson95"][0]>=frozen["positive_gate_Wilson95_lower_min"])
    return truths,passed


def run(args):
    contract=json.loads(CONTRACT.read_text())
    if contract["status"]!="FROZEN_BEFORE_FIRST_NUMERIC_RATE_OPENING":
        raise ValueError("First-shot protocol not frozen")
    original=json.loads(NATIVE_LOCK.read_text())["native_original_transport"]
    assets=json.loads(CLIMATE.read_text())["assets"]
    args.work.mkdir(parents=True,exist_ok=True);args.output.parent.mkdir(parents=True,exist_ok=True)
    out={"schema":"ttf_PROCEED70_Haldane_source_only_known_truth_calibration_v01",
         "actual_genetic_rate_mean_SD_or_fitness_outcome_parsed":False,
         "real_Haldane_or_Darwin_outcome_fitted":False,
         "prior_TTF_frozen_outcomes_reopened":False,
         "raw_study_species_site_or_biological_response_exported":False}
    try:
        catalog=validated_get_catalog("https://borealisdata.ca/api/datasets/:persistentId/?persistentId=doi%3A10.5683%2FSP3%2FNXSL3Q")
        verify_frozen_catalog(catalog,json.loads(SOURCE_CONTRACT.read_text())["current_source"])
        raw,md5,sha=fetch_original_datafile_bytes(1024206)
        if sha!=original["observed_sha256"] or md5!=original["observed_md5"] or len(raw)!=original["observed_size_bytes"]:
            raise ValueError("Source-only native original PROCEED MD5/SHA drift")
        full,studies=source_geometry(raw,json.loads(FULL_GEO.read_text()))
        coords=sorted({v for _,_,a,b in full for v in (a,b)})
        xy=np.asarray(coords,float)
        if hashlib.sha256(xy.astype("<f8").tobytes()).hexdigest()!="aed6486ee9aa49b0d454426c4d5964ebd79a8a96db2aac452c80fd5b7f52af32":
            raise ValueError("Original 253-source-coordinate SHA mismatch")
        climate=[];receipts=[]
        for i,a in enumerate(assets):
            file=args.work/f"source_{i}.tif"
            verified_raster(a,file)
            vals=raster_values(file,xy)
            climate.append(vals);file.unlink()
            receipts.append({"bioclim":f"{a['variable']}_{a['time_index']}",
                             "SHA256":a["sha256"],"finite_original_sites":int(np.isfinite(vals).sum())})
            print("HALDANE_SOURCE_ONLY_CHELSA_VERIFIED",i+1,flush=True)
        matrix=np.column_stack(climate)
        if matrix.shape!=(253,8) or not np.isfinite(matrix).all() or (
            hashlib.sha256(matrix.astype("<f8").tobytes()).hexdigest()
             !="7d03b2f158ade223b752a2b7d6ac57c38fd163e0b44f5566a01e8e6fe6d6b4fb"):
            raise ValueError("Focal Haldane source-only climate raster matrix unreproducible")
        X,H,sp,st,role,climatesd,mask=source_only_geometry(raw,matrix,coords,json.loads(FULL_GEO.read_text()))
        info=source_history_concentration(H,X,sp,role)
        c=contract["synthetic_pre_response_qualification"]
        qualified=(info["history_unique_fraction_after_all_21_controls"]>=c["source_train_history_unique_fraction_min"]
                   and info["effective_source_species_of_history_residual"]>=c["source_train_history_effective_species_min"]
                   and info["maximum_one_source_species_history_residual_share"]<=c["max_one_source_species_fraction_of_history_residual_information"])
        out.update(info)
        out.update({"source_8_chelsa_sha":"7d03b2f158ade223b752a2b7d6ac57c38fd163e0b44f5566a01e8e6fe6d6b4fb",
                    "source_original_sites":253,"source_Haldane_presence_unique_geo_pairs":len(H),
                    "source_Haldane_original_species_train":29,"source_Haldane_original_species_heldout":27,
                    "source_train_0BP_climate_SD":climatesd.tolist(),
                    "source_model_21controls_full_numeric_sha256":hashlib.sha256(X.astype("<f8").tobytes()).hexdigest(),
                    "source_history_magnitude_sha256":hashlib.sha256(H.astype("<f8").tobytes()).hexdigest(),
                    "source_8_assets":receipts})
        if not qualified:
            out["status"]="NOT_EVALUABLE_HALDANE_SOURCE_CONCENTRATION_GATE"
        else:
            null,passed=known_truth_gate(X,H,sp,st,role,contract)
            out["known_truth"]=null
            out["status"]="PASS_HALDANE_SYNTHETIC_KNOWN_TRUTH_BEFORE_RATE_OPENING" if passed else c["if_fail"]
    except Exception as exc:
        out["status"]="NOT_EVALUABLE_PRE_RESPONSE_SOURCE_OR_SYNTHETIC_GATE"
        out["error_type"]=type(exc).__name__
        out["error_message"]=str(exc)[:400]
    args.output.write_text(json.dumps(out,sort_keys=True,indent=2)+"\n")
    print(json.dumps({k:out.get(k) for k in (
        "status","history_unique_fraction_after_all_21_controls",
        "effective_source_species_of_history_residual",
        "maximum_one_source_species_history_residual_share",
        "source_Haldane_presence_unique_geo_pairs")},sort_keys=True))


if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--work",type=Path,required=True)
    p.add_argument("--output",type=Path,required=True)
    run(p.parse_args())
