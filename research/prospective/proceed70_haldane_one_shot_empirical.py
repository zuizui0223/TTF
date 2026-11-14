"""One authorized prospective empirical Haldane-rate prediction, source-gated.

DO NOT INVOKE until synthetic source-only and exact original CHELSA8 audits PASS.
This is a new observational source category of genetically influenced
phenotypic change rates; it does NOT identify Darwinian local adaptation,
reaction-norm plasticity or species diversification.

Licensed/source-local individual rates and coordinates stay in the runner;
only aggregate study/species statistics and source hashes can be uploaded.
"""
from __future__ import annotations
import argparse,csv,hashlib,io,json,math
from collections import defaultdict,Counter
from pathlib import Path

import numpy as np
from proceed_v62_metadata_row_census import (
    fetch_original_datafile_bytes,valid_binomial,valid_coordinates)
from proceed_v62_borealis_header_only import validated_get_catalog,verify_frozen_catalog
from proceed_v62_genetic_geometry_intersection import SOURCE_CONTRACT,NATIVE_LOCK
from proceed70_haldane_presence_mask import (CONTRACT as PRESENCE,FIELDS,META,
                                             is_present_raw_token,missingness_summary)
from proceed70_haldane_source_null_calibration import (
    CONTRACT as FULLCONTRACT,CLIMATE,FULL_GEO,presence_keyed_metadata,
    source_only_geometry
)
from proceed_plant_chelsa_source_support import source_geometry,verified_raster
from proceed70_haldane_rate_model_core import (
    ridge_predict,species_mse_comparison,SEED
)
from site_feature_core import raster_values

HERE=Path(__file__).resolve().parent


def author_qualified_Haldane(row):
    """Original v6.2 BasicCode lines 85-124, plus explicitly frozen finite guards.

    Numeric trait means, SDs, sizes are read ONLY in this outcome-opening
    function. All code/source/predictor/synthetic gates must pass first.
    """
    try:
        a,b,sd1,sd2,n1,n2,g=[float(row[k]) for k in FIELDS]
    except (KeyError,ValueError,TypeError):
        return None
    vals=(a,b,sd1,sd2,n1,n2,g)
    if not all(math.isfinite(z) for z in vals) or min(sd1,sd2)<=0 or min(n1,n2)<2 or g<=0:
        return None
    # Frozen official source outlier CV bounds; zero means imply an infinite CV
    # and fail rather than inventing an epsilon or silently changing endpoint.
    if a==0 or b==0:return None
    cv1=100*abs(sd1/a);cv2=100*abs(sd2/b)
    if not(0.1<=cv1<=1000 and 0.1<=cv2<=1000):return None
    dof=n1+n2-2
    if dof<=0:return None
    pooled2=((n1-1)*sd1**2+(n2-1)*sd2**2)/dof
    if not math.isfinite(pooled2) or pooled2<=0:return None
    H=(b-a)/math.sqrt(pooled2)/g
    if not math.isfinite(H):return None
    return float(math.log1p(abs(H)))


def only_after_authorized_response(raw,source_keys,authority):
    """Returns private aggregated outcomes keyed by fixed pre-response source geography."""
    if authority.get("status")!="PASS_HALDANE_SYNTHETIC_KNOWN_TRUTH_BEFORE_RATE_OPENING":
        raise ValueError("No authorized single biological response opening")
    p=json.loads(PRESENCE.read_text())
    missing={x.lower() for x in p["field_presence_only"]["absent_if_trimmed_lowercase_in"]}
    r=csv.reader(io.TextIOWrapper(io.BytesIO(raw),encoding="utf-8-sig",newline=""),
                 delimiter=",",strict=True)
    names=next(r)
    columns=list(META)+list(FIELDS)
    if len(names)!=69 or len(set(names))!=69 or not set(columns).issubset(names):
        raise ValueError("Source before response opening not original 69-column data")
    idx={key:names.index(key) for key in columns}
    all_numeric_source_rows=0
    n_mask_rows=0
    count_valid=0
    private=defaultdict(list)
    for values in r:
        if not values:continue
        if len(values)!=69:raise ValueError("Source width drift after authorization")
        m={k:values[idx[k]].strip() for k in META}
        if not(m["phylum"]=="Streptophyta" and m["taxa"]=="Plant" and m["genphen"]=="Genetic"):
            continue
        a=valid_coordinates(m["sample1_latitude"],m["sample1_longitude"])
        b=valid_coordinates(m["sample2_latitude"],m["sample2_longitude"])
        species=valid_binomial(m["sp_pub"]) or valid_binomial(m["sp_ncbi"])
        study=m["ref_ID"]
        if not(a and b and species and study):continue
        all_numeric_source_rows+=1
        if not all(is_present_raw_token(values[idx[x]],missing) for x in FIELDS):
            continue
        n_mask_rows+=1
        key=(species,study,(a[1],a[0]),(b[1],b[0]))
        if key not in source_keys:
            raise ValueError("Outcome opening source geometry differs from frozen Haldane mask")
        # FIRST permitted numeric rate-field access happens here.
        rate=author_qualified_Haldane({k:values[idx[k]] for k in FIELDS})
        if rate is None:continue
        private[key].append(rate)
        count_valid+=1
    if all_numeric_source_rows!=1561 or n_mask_rows!=1279:
        raise ValueError("Response input source panel no longer frozen")
    return private,{"source_plant_eligible_records":all_numeric_source_rows,
                    "H_presence_mask_rows":n_mask_rows,
                    "canonical_original_H_numeric_valid_rows":count_valid,
                    "canonical_H_rejected_after_presence_mask":n_mask_rows-count_valid}


def heldout_first_shot(raw,source_info,contract):
    meta=presence_keyed_metadata(raw,json.loads(PRESENCE.read_text()))
    keys=sorted(meta)
    if len(keys)!=165:
        raise ValueError("H source geometry changed before response opening")
    # Numeric input authority checks must already have passed.
    private,qc=only_after_authorized_response(raw,set(keys),source_info)
    allowed=sorted(private)
    roles=source_info["_private_roles"]
    X,H=source_info["_private_X"],source_info["_private_H"]
    chosen=np.array([key in private for key in keys])
    sp=np.array([key[0] for key in keys]);study=np.array([key[1] for key in keys])
    role=np.array([roles[s] for s in sp])
    y=np.asarray([float(np.median(private[key])) for key in keys if key in private])
    used_X=X[chosen];used_H=H[chosen];used_sp=sp[chosen];used_study=study[chosen];used_role=role[chosen]
    bounds=contract["biology"]["final_first_response_QC_gate"]
    source_species=set(used_sp)
    train=set(used_sp[used_role=="train"]);heldout=set(used_sp[used_role=="heldout"])
    if (len(source_species)<bounds["min_taxa"] or len(train)<bounds["min_train_taxa"]
        or len(heldout)<bounds["min_heldout_taxa"] or
        len(set(used_study))<bounds["min_source_studies"] or
        len({s.split()[0] for s in source_species})<bounds["min_genera"] or
        int((used_role=="train").sum())<bounds["min_train_pair_geometries"] or
        int((used_role=="heldout").sum())<bounds["min_heldout_pair_geometries"]):
        return {
            "status":bounds["fail"],
            "empirical_source_rate_numeric_QC_opened":True,
            "empirical_historical_predictive_effect_tested":False,
            "numeric_Haldane_QC":qc,
            "final_species_count":len(source_species),
            "final_train_species":len(train),"final_heldout_species":len(heldout),
            "final_original_studies":len(set(used_study)),
            "final_unique_geographic_comparison_pairs":len(y)
        }
    b=ridge_predict(used_X,y,used_sp,used_role)
    h=ridge_predict(used_X,y,used_sp,used_role,history=used_H)
    res=species_mse_comparison(y,b,h,used_sp,used_role,reps=1999,seed=SEED)
    return {
        "status":"QUALIFIED_POSITIVE_HISTORICAL_PREDICTION_OF_GENETICALLY_BASED_TRAIT_CHANGE_RATE" if res["positive"] else
                  "NO_QUALIFIED_ADDED_HISTORICAL_PREDICTION_OF_TRAIT_CHANGE_RATE",
        "empirical_source_rate_numeric_QC_opened":True,
        "empirical_historical_predictive_effect_tested":True,
        "numeric_Haldane_QC":qc,
        "final_species_count":len(source_species),
        "final_train_species":len(train),"final_heldout_species":len(heldout),
        "final_original_studies":len(set(used_study)),
        "final_unique_geographic_comparison_pairs":len(y),
        "source_study_species_joint_holdout_preserved":not bool(set(used_study[used_role=="train"])&set(used_study[used_role=="heldout"])),
        "species_equal_heldout_evaluation":res,
        "gene_by_environment_local_adaptation_causality_identified":False,
        "plasticity_reaction_norm_identified":False,
        "speciation_reproductive_isolation_identified":False
    }


def run(args):
    contract=json.loads(FULLCONTRACT.read_text())
    if contract["status"]!="FROZEN_BEFORE_FIRST_NUMERIC_RATE_OPENING":
        raise ValueError("Biological outcome contract not frozen")
    handoff=json.loads(args.handoff.read_text())
    if handoff.get("status")!="PASS_HALDANE_SYNTHETIC_KNOWN_TRUTH_BEFORE_RATE_OPENING":
        raise ValueError("Pre-response synthetic gate not authorized")
    proof=contract["synthetic_pre_response_qualification"]
    if not(handoff["known_truth"]["null"]["Wilson95"][1]<=proof["null_gate_Wilson95_upper_max"]
        and handoff["known_truth"]["positive"]["Wilson95"][0]>=proof["positive_gate_Wilson95_lower_min"]):
        raise ValueError("Synthetic null/power gate disagrees with authority")
    native=json.loads(NATIVE_LOCK.read_text())["native_original_transport"]
    src=json.loads(SOURCE_CONTRACT.read_text())["current_source"]
    args.work.mkdir(parents=True,exist_ok=True);args.output.parent.mkdir(parents=True,exist_ok=True)
    output={"schema":"ttf_proceed70_first_shot_genetic_based_plant_phenotypic_Haldane_rate_history_prediction_v01",
            "status":"SOURCE_OR_KNOWN_TRUTH_AUTHORITY_STOP_NO_RATE_OPENING",
            "empirical_source_rate_numeric_QC_opened":False,
            "empirical_historical_predictive_effect_tested":False,
            "prior_TTF_Q_COI_He_Fst_Jensen_outcomes_unchanged":True,
            "raw_study_site_climate_or_genetic_rate_values_exported":False}
    try:
        catalog=validated_get_catalog("https://borealisdata.ca/api/datasets/:persistentId/?persistentId=doi%3A10.5683%2FSP3%2FNXSL3Q")
        verify_frozen_catalog(catalog,src)
        raw,md5,sha=fetch_original_datafile_bytes(1024206)
        if (md5!=native["observed_md5"] or sha!=native["observed_sha256"] or
            len(raw)!=native["observed_size_bytes"]):
            raise ValueError("Source original MD5/SHA identity changed")
        full,_=source_geometry(raw,json.loads(FULL_GEO.read_text()))
        coords=sorted({v for _,_,a,b in full for v in (a,b)})
        xy=np.asarray(coords,float)
        if hashlib.sha256(xy.astype("<f8").tobytes()).hexdigest()!=contract["evidence_chain"]["original_sites_SHA256"]:
            raise ValueError("Original source coordinate SHA not reproduced")
        assets=json.loads(CLIMATE.read_text())["assets"]
        source_vals=[]
        for i,a in enumerate(assets):
            p=args.work/f"original_{i}.tif"
            verified_raster(a,p)
            source_vals.append(raster_values(p,xy))
            p.unlink()
            print("HALDANE_FIRST_SHOT_PREOPEN_CLIMATE_SOURCE",i+1,flush=True)
        x=np.column_stack(source_vals)
        if hashlib.sha256(x.astype("<f8").tobytes()).hexdigest()!=contract["evidence_chain"]["original_chelsa8_SHA256"]:
            raise ValueError("Source CHELSA numeric original SHA not reproducible")
        X,H,sp,st,role,sd,mask=source_only_geometry(raw,x,coords,json.loads(FULL_GEO.read_text()))
        if (hashlib.sha256(X.astype("<f8").tobytes()).hexdigest()!=handoff["source_model_21controls_full_numeric_sha256"] or
            hashlib.sha256(H.astype("<f8").tobytes()).hexdigest()!=handoff["source_history_magnitude_sha256"]):
            raise ValueError("Source model/history numeric pre-response SHA does not match frozen known-truth")
        # FIREWALL: only AFTER exact original PROCEED, CHELSA, modeled source
        # and frozen known-truth qualification can original mean, SD and n
        # magnitudes be used in the next single function.
        print("AUTHORIZED_FIRST_ONLY_PROCEED_H_RATE_RESPONSE_OPENING",flush=True)
        rs=source_geometry(raw,json.loads(FULL_GEO.read_text()))[1]
        from proceed_plant_chelsa_source_support import component_assignment
        roles=component_assignment(rs)[0]
        private={"_private_roles":roles,"_private_X":X,"_private_H":H,
                 "status":handoff["status"]}
        output.update(heldout_first_shot(raw,private,contract))
        output.update({"input_21_baseline_SHA256":handoff["source_model_21controls_full_numeric_sha256"],
                       "input_single_history_SHA256":handoff["source_history_magnitude_sha256"],
                       "original_rate_source_SHA256":sha,
                       "original_climate_matrix_SHA256":contract["evidence_chain"]["original_chelsa8_SHA256"],
                       "source_only_synthetic_gate_run":handoff["source_only_run_id"],
                       "frozen_first_shot_contract":"research/prospective/proceed70_haldane_first_shot_contract_v01.json"})
    except Exception as exc:
        if output.get("empirical_source_rate_numeric_QC_opened"):
            output["status"]="TERMINAL_EMPIRICAL_OUTCOME_OPENED_ERROR_NO_RERUN_OR_RETUNING"
        else:
            output["status"]="SOURCE_OR_KNOWN_TRUTH_AUTHORITY_STOP_NO_RATE_OPENING"
        output["error_type"]=type(exc).__name__
        output["error_message"]=str(exc)[:400]
    args.output.write_text(json.dumps(output,sort_keys=True,indent=2)+"\n")
    print(json.dumps({k:output.get(k) for k in (
        "status","empirical_source_rate_numeric_QC_opened",
        "empirical_historical_predictive_effect_tested",
        "final_species_count","final_train_species","final_heldout_species",
        "species_equal_heldout_evaluation")},sort_keys=True))


if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--handoff",required=True,type=Path)
    p.add_argument("--work",required=True,type=Path)
    p.add_argument("--output",required=True,type=Path)
    run(p.parse_args())
