"""Run Stage21 exact-site synthetic spatial-null/power calibration. No empirical AUC."""
import argparse,hashlib,json
from pathlib import Path
import numpy as np
from spatial_null_core import (BASELINE_COEFF,standardize_predictors,
    training_history_unique_fraction,draw_synthetic_labels,evaluate_species_holdout)
EXPECTED={"coordinates_raw_sha256":"fad49f6ec345ef45c1387ef76c8b3b724caf5f39322f8e9c3da56227a7a40bfd",
 "feature_matrix_raw_sha256":"9c408f6f23be75bee6aa79475bc68542ade225927079e7a43548dc152cdc2db8",
 "class_labels_sha256":"e08632a14a990804a7a1c83e236b70834056ac61b7157b9832b8f91c0df855d5"}
NAMES=["bio01_-190","bio07_-190","bio12_-190","bio15_-190",
       "bio01_20","bio07_20","bio12_20","bio15_20",
       "mountain_km","water_km","ecoregion1_km","ecoregion2_km","TRI"]
def wilson95(k,n):
    z=1.959963984540054;p=k/n;d=1+z*z/n
    c=(p+z*z/(2*n))/d;h=z*np.sqrt(p*(1-p)/n+z*z/(4*n*n))/d
    return [float(max(0.,c-h)),float(min(1.,c+h))]

def check_receipt_hashes(receipt,expected):
    if any(receipt.get(k)!=v for k,v in expected.items()):
        raise ValueError("replayed source receipt does not match chosen frozen numerical input")

def locked_site_arrays(private,receipt_path,panel_path,expected):
    r=json.loads(Path(receipt_path).read_text())
    if r["status"]!="PASS_INPUT_SUPPORT_NOT_EMPIRICAL_EFFECT" or r["analysis_boundary"]["genetic_model_fitted"]:
        raise ValueError("Stage20 input firewall violated")
    check_receipt_hashes(r,expected)
    with np.load(private,allow_pickle=False) as z:
        xy=np.asarray(z["xy"],float);x=np.asarray(z["features"],float)
        case=np.asarray(z["case"],dtype="u1");species=np.asarray(z["species_index"],dtype=int)
        names=z["feature_names"].tolist()
    if names!=NAMES or len(x)!=76700 or len(xy)!=76700 or case.sum()!=10013:
        raise ValueError("frozen features/labels changed")
    for array,k in ((xy.astype("<f8"),"coordinates_raw_sha256"),
                    (x.astype("<f8"),"feature_matrix_raw_sha256"),(case,"class_labels_sha256")):
        if hashlib.sha256(array.tobytes()).hexdigest()!=expected[k]:
            raise ValueError("actual private feature matrix SHA does not match original")
    panel=json.loads(Path(panel_path).read_text())
    kept=[s for s in panel["qualified_species_records"] if s["synthetic_estimator_min10_cases"]]
    if len(kept)!=80 or not np.array_equal(np.unique(species),np.arange(80)):
        raise ValueError("wrong species identification")
    roles_by_id=np.asarray([s["split_role"] for s in kept],dtype=int)
    if np.bincount(roles_by_id,minlength=2).tolist()!=[44,36]:
        raise ValueError("heldout roles changed")
    roles=roles_by_id[species]
    counts={i:int(case[species==i].sum()) for i in range(80)}
    if sum(counts.values())!=10013 or min(counts.values())<10:
        raise ValueError("frozen source case-support changed")
    return x,xy,species,roles,counts

def run(args):
    c=json.loads(args.contract.read_text())
    if c["status"]!="FROZEN_BEFORE_ANY_GENETIC_HISTORY_AUC":
        raise ValueError("spatial-null protocol not locked")
    if not np.array_equal(BASELINE_COEFF,np.asarray(c["baseline_coefficients_in_training_SD_units"])):
        raise ValueError("synthetic baseline coefficient drift")
    amend=json.loads(args.amendment.read_text())
    if amend.get("status")!="TWO_INDEPENDENT_REPLAYS_MATCH_PRE_BIOLOGICAL_AUC":
        raise ValueError("No authorized numerical source replay amendment")
    if amend["original_run"]["matrix_sha256"]!=EXPECTED["feature_matrix_raw_sha256"]:
        raise ValueError("Missing original numerical SHA discrepancy provenance")
    new_sha=amend["new_frozen_feature_matrix_sha256"]
    if amend["replay_one"]["feature_matrix_sha256"]!=new_sha or amend["replay_two"]["feature_matrix_sha256"]!=new_sha:
        raise ValueError("Two independent original-source replay SHAs disagree")
    expected=dict(EXPECTED)
    expected["feature_matrix_raw_sha256"]=new_sha
    r1=json.loads(args.reference_stage21.read_text())
    r2=json.loads(args.reference_stage20.read_text())
    for ref in (r1,r2):
        for key,val in expected.items():
            if ref.get(key)!=val:raise ValueError("Replay receipt disagrees with amended source")
        if ref.get("status")!="PASS_INPUT_SUPPORT_NOT_EMPIRICAL_EFFECT":
            raise ValueError("A replay failed its own original source gate")
    for key in ("source_receipts","species_support","complete_counts_by_feature","software"):
        if r1.get(key)!=r2.get(key):raise ValueError("Independent source receipts disagree: "+key)
    x,xy,sp,role,counts=locked_site_arrays(args.private,args.stage20_receipt,args.panel,expected)
    p,mu,scale,climate_sd=standardize_predictors(x,xy,sp,role)
    unique=training_history_unique_fraction(p,sp,role)
    required=c["synthetic_gate"]["min_train_historical_unique_fraction"]
    result={"schema":"ttf_jensen80_exact_design_spatial_null_known_truth_v01",
            "training_history_unique_fraction":unique,
            "frozen_unique_minimum":required,
            "train_species":44,"heldout_species":36,
            "observed_case_counts_used_only_as_simulation_fixed_margins":True,
            "case_site_records":10013,"background_site_records":66687,
            "genetic_history_auc_computed":False,"other_genetic_response_opened":False,
            "input_coordinate_sha256":EXPECTED["coordinates_raw_sha256"],
            "input_feature_sha256":expected["feature_matrix_raw_sha256"],
            "input_original_stage20_unreproduced_feature_sha256":EXPECTED["feature_matrix_raw_sha256"],
            "input_replay_amendment_status":amend["status"],
            "input_case_label_sha256":EXPECTED["class_labels_sha256"],
            "historical_training_current_bioclim_sd":climate_sd.tolist(),
            "positional_uncertainty_50km_gate":"SEPARATELY_REQUIRED"}
    if unique<required:
        result["status"]="NOT_EVALUABLE_HISTORICAL_INFORMATION_GATE"
    else:
        cfg=c["synthetic_generating_process"];rep=c["bootstrap"]["draws_per_world"]
        for kind,coefficient in [("null",cfg["null_history_log_odds"]),
                                 ("positive",cfg["positive_history_log_odds_per_training_SD"])]:
            trials=[]
            for i in range(cfg["reps_"+kind]):
                seed=cfg["seed"]+i+(100000 if kind=="positive" else 0)
                y=draw_synthetic_labels(p,xy,sp,counts,np.random.default_rng(seed),
                    history_log_odds=coefficient,spatial_sigma=cfg["spatial_random_field_amplitude"])
                result0=evaluate_species_holdout(p,y,sp,role,bootstrap_reps=rep,seed=seed+500000)
                trials.append({"index":i,"delta_auc":result0["mean_delta_auc"],
                    "ci_9833":result0["ci_9833"],
                    "positive":result0["positive_decision"]})
                if (i+1)%25==0:print("SYNTHETIC_PROGRESS",kind,i+1,flush=True)
            k=sum(t["positive"] for t in trials);n=len(trials)
            result[kind]={"worlds":n,"positive_decisions":k,
              "rate":k/n,"wilson95":wilson95(k,n),
              "median_delta_auc":float(np.median([t["delta_auc"] for t in trials])),
              "worlds_aggregated":trials}
        g=c["synthetic_gate"]
        passed=(result["null"]["wilson95"][1]<=g["max_null_Wilson95_upper"] and
                result["positive"]["wilson95"][0]>=g["min_positive_Wilson95_lower"])
        result["status"]="PASS_SYNTHETIC_SPATIAL_GATE" if passed else "NOT_EVALUABLE_SYNTHETIC_SPATIAL_GATE"
        result["limitations"]="Spatial Fourier Gumbel case-count worlds are conditional simulations, not reconstructions of genetic break line generation; full 50km positional sensitivity still required."
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(result,sort_keys=True,indent=2)+"\n")
    print(json.dumps({k:v for k,v in result.items() if k not in ("null","positive")},sort_keys=True))
if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--private",type=Path,required=True)
    p.add_argument("--stage20-receipt",type=Path,required=True)
    p.add_argument("--panel",type=Path,required=True)
    p.add_argument("--contract",type=Path,required=True)
    p.add_argument("--amendment",type=Path,required=True)
    p.add_argument("--reference-stage21",type=Path,required=True)
    p.add_argument("--reference-stage20",type=Path,required=True)
    p.add_argument("--output",type=Path,required=True)
    run(p.parse_args())
