"""Train-only climate-predictor identifiability; independent of Fst values."""
import argparse,csv,hashlib,json,math
from pathlib import Path
import numpy as np

def km(a,b):
    a0,a1,b0,b1=map(math.radians,(*a,*b))
    u=math.sin((b0-a0)/2)**2+math.cos(a0)*math.cos(b0)*math.sin((b1-a1)/2)**2
    return 12742.0176*math.asin(min(1.,math.sqrt(u)))
def read_climate(root):
    report=json.loads((root/"climate_receipt.json").read_text())
    pairs=json.loads((root/"eligible_pair_coordinates.json").read_text())
    if report["status"]!="PASS_CLIMATE_RESPONSE_BLIND" or pairs["status"]!=report["status"]:
        raise ValueError("No climate authorization")
    if report["Fst_values_opened"] or pairs["Fst_values_opened"]:
        raise ValueError("Fst response firewall already open")
    climate=np.load(root/"climate_values.npy",allow_pickle=False)
    if hashlib.sha256(climate.tobytes()).hexdigest()!=report["climate_value_matrix_sha256"]:
        raise ValueError("climate bytes mismatch")
    with (root/"site_keys.csv").open(newline="") as f:
        z=list(csv.DictReader(f))
    sites=[(float(s["latitude"]),float(s["longitude"])) for s in z]
    if len(sites)!=len(climate) or len(set(sites))!=len(sites):
        raise ValueError("site-index mismatch")
    return report,pairs["groups"],dict(zip(sites,climate))

def features(group,by_site,sd):
    vals=[]
    for lat1,lon1,lat2,lon2 in group["climate_valid_pair_coordinates"]:
        a=(lat1,lon1);b=(lat2,lon2)
        v=by_site[a];w=by_site[b]
        if not np.isfinite(v).all() or not np.isfinite(w).all():
            raise ValueError("climate-missing pair in eligible list")
        cur=float(np.linalg.norm((v[4:]-w[4:])/sd))
        dh=(v[:4]-v[4:])-(w[:4]-w[4:])
        hist=float(np.linalg.norm(dh/sd))
        vals.append([np.log1p(km(a,b)),cur,(a[0]+b[0])/2,
                    (a[1]+b[1])/2,hist])
    x=np.asarray(vals,float)
    if x.ndim!=2 or x.shape[0]<5 or x.shape[1]!=5 or not np.isfinite(x).all():
        raise ValueError("invalid study-group covariates")
    return x-x.mean(axis=0)

def partial_fraction(groups):
    X=[];y=[]
    for x in groups:
        w=1/math.sqrt(len(x))
        X.append(w*x[:,:4]);y.append(w*x[:,4])
    a=np.vstack(X);b=np.concatenate(y)
    norm=float(b@b)
    if norm<=1e-12:raise ValueError("no historical predictor signal")
    coef=np.linalg.lstsq(a,b,rcond=1e-12)[0]
    residue=b-a@coef
    return float((residue@residue)/norm)

def gate(root):
    report,groups,by_site=read_climate(root)
    tr={site for g in groups if g["role"]==0
        for p in g["climate_valid_pair_coordinates"]
        for site in ((p[0],p[1]),(p[2],p[3]))}
    pts=np.asarray([by_site[s][4:] for s in tr])
    sd=pts.std(axis=0)
    if np.any(sd<=1e-10):raise ValueError("training climate variable constant")
    train=[features(g,by_site,sd) for g in groups if g["role"]==0]
    frac=partial_fraction(train)
    ok=frac>=.05
    return {"schema":"ttf_calipopgen_fst_train_climate_information_v01",
      "status":"PASS_PREDICTOR_INFORMATION_BEFORE_FST" if ok else "NOT_EVALUABLE_HISTORY_REDUNDANT",
      "train_climate_sites":len(tr),"training_group_count":len(train),
      "unique_historical_predictor_fraction":frac,
      "frozen_minimum_fraction":.05,
      "climate_input_sha256":report["climate_value_matrix_sha256"],
      "baseline":["log1p_geodesic_distance","0BP_multivariable_distance","midpoint_lat","midpoint_lon"],
      "focal":"21ka-minus-0BP_displacement_vector_difference_norm",
      "training_current_bioclim_sd":sd.tolist(),
      "Fst_numeric_values_opened":False}
if __name__=="__main__":
    a=argparse.ArgumentParser();a.add_argument("--climate-dir",type=Path,required=True)
    a.add_argument("--output",type=Path,required=True)
    v=a.parse_args()
    result=gate(v.climate_dir)
    v.output.parent.mkdir(parents=True,exist_ok=True)
    v.output.write_text(json.dumps(result,indent=2,sort_keys=True)+"\n")
    print(json.dumps(result,sort_keys=True))
