"""One-shot biological opening for GenDivRange; run only after frozen climate gate."""
import argparse,collections,csv,hashlib,json,math
from pathlib import Path
import numpy as np
from inspect_headers_only import fetch_verified
from reconcile_keys import stream
from model_he_holdout import evaluate

def run(args):
    auth=json.loads(args.authorization.read_text())
    if auth.get("status")!="AUTHORIZE_NEW_GENDIVRANGE_HE_ONLY":
        raise ValueError("no preregistered He response authorization")
    root=args.climate_dir
    geo=json.loads((root/"eligible_geometry.json").read_text())
    receipt=json.loads((root/"results/gendivrange-climate/climate_receipt.json").read_text())
    matrix=np.load(root/"results/gendivrange-climate/climate_values.npy",allow_pickle=False)
    if hashlib.sha256(matrix.tobytes()).hexdigest()!=auth["climate_matrix_raw_sha256"]:
        raise ValueError("frozen exact-CHELSA climate matrix mismatch")
    if receipt["genetic_He_values_opened"] or matrix.shape[1]!=8:
        raise ValueError("unexpected climate provenance")
    with (root/"results/gendivrange-climate/site_keys.csv").open(newline="") as f:
        sites=list(csv.DictReader(f))
    assert len(sites)==len(matrix)==auth["total_climate_sites"]
    manifest={(m["species_id"],m["study_id"]):m for m in geo["adult_terrestrial_frozen_manifest"]}
    lookup={}
    for r,v in zip(sites,matrix):
        key=(r["species_id"],r["study_id"],float(r["latitude"]),float(r["longitude"]))
        if key in lookup:raise ValueError("duplicate frozen climate site")
        if np.isfinite(v).all():lookup[key]=v
    # Genetics are opened only below this exact fence.
    meta=json.loads(args.metadata.read_text())
    fs={x["filename"]:x for x in meta["files"]}
    pop=fetch_verified(fs["pop_tab_v2025-03-31.xlsx"])
    group_pop=collections.defaultdict(lambda:collections.defaultdict(list))
    excluded=collections.Counter()
    for r in stream(pop,{"Spec_id","Study_id","Pop_id","Latitude","Longitude","He","N"}):
        k=(r.get("Spec_id"),r.get("Study_id"))
        if k not in manifest:continue
        try:
            lat=float(r["Latitude"]);lon=float(r["Longitude"]);v=lookup.get((*k,lat,lon))
            if v is None:continue
            value=float(r["He"])
            if not 0<=value<=1:raise ValueError("He out of range")
        except (TypeError,ValueError,KeyError):
            excluded["invalid_or_missing_He"]+=1;continue
        popid=r.get("Pop_id")
        if not popid:excluded["no_population_id"]+=1;continue
        try:
            N=float(r.get("N"))
            if not np.isfinite(N) or N<=0:N=float("nan")
        except (ValueError,TypeError):N=float("nan")
        group_pop[k][str(popid)].append((lat,lon,value,N))
    groups={}
    for k,byid in group_pop.items():
        bysite=collections.defaultdict(list)
        for pid,items in byid.items():
            values={p[2] for p in items}
            positions={(p[0],p[1]) for p in items}
            if len(values)!=1 or len(positions)!=1:
                excluded["conflicting_duplicate_population_id"]+=1;continue
            if len(items)>1:excluded["identical_duplicate_rows"]+=len(items)-1
            x=items[0]
            bysite[(x[0],x[1])].append((x[2],x[3]))
        X=[];Y=[]
        for (lat,lon),recs in sorted(bysite.items()):
            climate=lookup[(*k,lat,lon)]
            H=float(np.mean([v[0] for v in recs]))
            ns=[v[1] for v in recs if np.isfinite(v[1]) and v[1]>0]
            N=float(np.exp(np.mean(np.log(ns)))) if ns else float("nan")
            X.append([*climate,lat,lon,N]);Y.append(H)
        if len(Y)>=5:
            m=manifest[k]
            groups[k]={"species":m["species_name"],"role":m["role"],
                       "x":np.asarray(X,float),"y":np.asarray(Y,float)}
    species={i:{g["species"] for g in groups.values() if g["role"]==i} for i in (0,1)}
    group_counts=[sum(g["role"]==i for g in groups.values()) for i in (0,1)]
    if len(species[0]|species[1])<40 or min(group_counts)<20 or species[0]&species[1]:
        raise RuntimeError("post-opening observed response support NOT_EVALUABLE")
    result=evaluate(groups,draws=1999,seed=20261029)
    result["observed_He_data_quality"]=dict(excluded)
    result["all_observed_He_species_groups"]=len(groups)
    result["source"]="GenDivRange v2025-03-31 Figshare 28407989"
    result["climate"]="SHA-256 pinned original C1 21ka and 0BP, eight CHELSA-TraCE21k COGs"
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(result,indent=2,sort_keys=True)+"\n")
    print(json.dumps(result,sort_keys=True))
if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--authorization",type=Path,required=True)
    p.add_argument("--metadata",type=Path,required=True)
    p.add_argument("--climate-dir",type=Path,required=True)
    p.add_argument("--output",type=Path,required=True)
    run(p.parse_args())
