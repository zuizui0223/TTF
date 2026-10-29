"""Exact-geometry synthetic null/power gate using shared-site pairwise dependence; NO Fst."""
import argparse,json,math
from pathlib import Path
import numpy as np
from climate_feature_gate import read_climate,features
from model_fst_holdout import evaluate

def wilson(k,n):
    z=1.959963984540054
    p=k/n
    d=1+z*z/n
    ctr=(p+z*z/(2*n))/d
    half=z*math.sqrt(p*(1-p)/n+z*z/(4*n*n))/d
    return [float(ctr-half),float(ctr+half)]

def get_geometry(root):
    _,groups,sites=read_climate(root)
    tr={pt for g in groups if g["role"]==0
        for p in g["climate_valid_pair_coordinates"]
        for pt in ((p[0],p[1]),(p[2],p[3]))}
    sd=np.std(np.asarray([sites[p][4:] for p in tr]),axis=0)
    X=[features(g,sites,sd) for g in groups]
    scale=np.sqrt(np.mean([np.mean(x*x,axis=0) for x,g in zip(X,groups)
                           if g["role"]==0],axis=0))
    if np.any(scale<=1e-9):raise ValueError("constant geometry feature")
    return groups,X,scale

def synthetic_world(groups,X,scale,settings,history_beta,rng):
    species=sorted({g["species"] for g in groups})
    private={s:rng.normal(0,settings["species_private_history_slope_sd"]) for s in species}
    baseline=np.asarray(settings["baseline_coefficients_in_training_column_standard_deviations"],float)
    out=[]
    for g,x in zip(groups,X):
        sites={tuple(p[:2]) for p in g["climate_valid_pair_coordinates"]}
        sites.update(tuple(p[2:]) for p in g["climate_valid_pair_coordinates"])
        latent={p:rng.normal(0,settings["population_latent_effect_sd"]) for p in sites}
        pair_effect=np.asarray([latent[tuple(p[:2])]+latent[tuple(p[2:])]
                       for p in g["climate_valid_pair_coordinates"]],float)
        b=history_beta+private[g["species"]]
        y=x[:,:4]/scale[:4]@baseline + b*x[:,4]/scale[4]+pair_effect
        y=y+rng.normal(0,settings["individual_dyad_noise_sd"],len(y))
        y-=y.mean()
        out.append({"species":g["species"],"role":g["role"],"x":x,"y":y})
    return out

def gate(root,contract):
    settings=contract["simulation"];q=contract["qualification"]
    groups,X,scale=get_geometry(root)
    samples={}
    for mode,n,beta in (
        ("null",settings["null_worlds"],settings["null_slope"]),
        ("positive",settings["positive_worlds"],settings["history_positive_slope_per_training_history_SD"])
    ):
        rejects=0;gains=[]
        for i in range(n):
            seed=int(settings["master_seed"])+i+(0 if mode=="null" else 100000)
            rng=np.random.default_rng(seed)
            synthetic=synthetic_world(groups,X,scale,settings,beta,rng)
            result=evaluate(synthetic,bootstrap_draws=settings["bootstraps_per_world"],seed=seed+500000)
            gains.append(result["relative_heldout_MSE_gain"])
            rejects+=int(result["positive_gain_95_lower_gt_zero"])
        samples[mode]={"worlds":n,"positive_decisions":rejects,
                       "positive_fraction":rejects/n,"wilson95":wilson(rejects,n),
                       "median_relative_mse_gain":float(np.median(gains))}
    good=(samples["null"]["wilson95"][1]<=q["null_max_wilson95_upper_falsepositive"]
         and samples["positive"]["wilson95"][0]>=q["positive_min_wilson95_lower_power"])
    return {"schema":"ttf_calipopgen_fst_known_truth_gate_result_v01",
            "status":"PASS_SYNTHETIC_BEFORE_FST" if good else "NOT_EVALUABLE_SYNTHETIC_GATE",
            "null":samples["null"],"positive":samples["positive"],
            "training_feature_scale":scale.tolist(),
            "input_climate_matrix_sha256":json.loads((root/"climate_receipt.json").read_text())["climate_value_matrix_sha256"],
            "Fst_numeric_values_opened":False,
            "no_retry_or_calibration_tuning_authorized":True}
if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--climate-dir",type=Path,required=True)
    p.add_argument("--contract",type=Path,required=True)
    p.add_argument("--output",type=Path,required=True)
    a=p.parse_args()
    model=json.loads(a.contract.read_text())
    if model["status"]!="LOCKED_BEFORE_FST_NUMERIC_RESPONSE":raise ValueError("unknown contract")
    out=gate(a.climate_dir,model)
    a.output.parent.mkdir(parents=True,exist_ok=True)
    a.output.write_text(json.dumps(out,indent=2,sort_keys=True)+"\n")
    print(json.dumps(out,sort_keys=True))
