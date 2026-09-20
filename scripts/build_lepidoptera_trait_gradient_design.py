#!/usr/bin/env python3
from __future__ import annotations
import argparse,csv,hashlib,json,zipfile,tempfile
from pathlib import Path
import numpy as np
from ttf.phylogatr_confirmatory import scan_phylogatr_phase1,choose_one_panel_per_species
from ttf.conditional_transfer import prepare_target_source_pools
from ttf.phylogatr_compact_execution import prepare_phylogatr_compact_ttf_design

ALIASES=("COI","CO1","COX1","COXI","CYTOCHROME C OXIDASE SUBUNIT I","CYTOCHROME C OXIDASE SUBUNIT 1")
DEFAULT_EXCLUSION=Path("docs/supporting/lepidoptera_prior_identity_opened_species_exclusion_v0.1.json")

def traits(path):
    out={}
    with path.open(newline="",encoding="utf-8-sig") as f:
        for x in csv.DictReader(f):
            n=x["Species"].strip()
            if n and n not in out: out[n]=x
    return out

def wing(x):
    keys=("WS_L","WS_U","FW_L","FW_U","WS_L_Fem","WS_U_Fem","WS_L_Mal","WS_U_Mal","FW_L_Fem","FW_U_Fem","FW_L_Mal","FW_U_Mal")
    vals=[float(x[k]) for k in keys if str(x.get(k,"")).strip() not in {"","NA"}]
    if not vals: raise ValueError("missing wing size")
    return float(np.mean(vals))
def habitat(x): return tuple(str(x[k]).strip() for k in ("CanopyAffinity","EdgeAffinity","MoistureAffinity","DisturbanceAffinity"))
def complete(x):
    try: wing(x); float(x["NumberOfHostplantFamilies"])
    except (ValueError,KeyError): return False
    return str(x.get("Voltinism","")).strip() not in {"","NA"} and all(v not in {"","NA"} for v in habitat(x))

def zscore(x):
    x=np.asarray(x,float); sd=np.std(x)
    return np.zeros_like(x) if sd<=np.sqrt(np.finfo(float).eps) else (x-np.mean(x))/sd

def trait_kernel(names,tr):
    w=np.asarray([wing(tr[n]) for n in names]); h=np.asarray([float(tr[n]["NumberOfHostplantFamilies"]) for n in names])
    wz=zscore(np.log1p(np.maximum(w,0))); hz=zscore(np.log1p(np.maximum(h,0)))
    kw=np.exp(-np.abs(wz[:,None]-wz[None,:])); kh=np.exp(-np.abs(hz[:,None]-hz[None,:]))
    vol=np.asarray([str(tr[n]["Voltinism"]).strip() for n in names]); kv=(vol[:,None]==vol[None,:]).astype(float)
    habitats=[habitat(tr[n]) for n in names]
    khab=np.mean(np.stack([(np.asarray([x[j] for x in habitats])[:,None]==np.asarray([x[j] for x in habitats])[None,:]).astype(float) for j in range(4)]),axis=0)
    k=(kw+kv+kh+khab)/4.0
    return k,{"wing_size":w.tolist(),"voltinism":vol.tolist(),"host_breadth":h.tolist(),"habitat":[list(x) for x in habitats]}

def geometry_kernel(names,geos):
    feat=[]
    for n in names:
        g=geos[n]; xyz=np.asarray(g.coordinates,float); centroid=np.mean(xyz,axis=0); extent=float(np.sqrt(np.sum(np.var(xyz,axis=0))))
        feat.append([*centroid,np.log1p(g.n_edges),np.log1p(g.n_localities),np.log1p(extent)])
    f=np.asarray(feat,float); z=np.column_stack([zscore(f[:,j]) for j in range(f.shape[1])])
    d2=np.mean((z[:,None,:]-z[None,:,:])**2,axis=2)
    return np.exp(-0.5*d2),f

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--archive",type=Path,required=True); ap.add_argument("--leptraits",type=Path,required=True)
    ap.add_argument("--exclusion-json",type=Path,default=DEFAULT_EXCLUSION); ap.add_argument("--output",type=Path,required=True)
    a=ap.parse_args()
    exclusion=json.loads(a.exclusion_json.read_text())
    if exclusion.get("schema")!="ttf_lepidoptera_prior_identity_opened_species_exclusion_v0.1": raise RuntimeError("exclusion schema drift")
    excluded=set(map(str,exclusion["species"]))|set(map(str,exclusion.get("additional_operator_exposure_exclusion",[])))
    with tempfile.TemporaryDirectory() as td:
        with zipfile.ZipFile(a.archive) as z:z.extractall(td)
        roots=[p.parent for p in Path(td).rglob("genes.txt") if (p.parent/"cite.txt").is_file()]
        if len(roots)!=1: raise RuntimeError("archive root drift")
        scan=scan_phylogatr_phase1(roots[0],aliases=ALIASES,excluded_species=excluded,min_localities=12,min_endpoint_training_edges=5,neighbor_fraction=.15)
        panels=choose_one_panel_per_species(scan.candidates)
    tr=traits(a.leptraits); chosen={p.species:p for p in panels if p.order=="Lepidoptera" and p.species in tr and complete(tr[p.species])}
    names=tuple(sorted(chosen,key=lambda n:(hashlib.sha256(f"lepidoptera-trait-transfer-v0.1|major_complete|{n}".encode()).hexdigest(),n)))
    if len(names)<180: raise RuntimeError("species-disjoint major-complete panel fell below predeclared feasibility floor")
    cut=len(names)//2; train=names[:cut]; evaluation=names[cut:]; geos={n:chosen[n].geometry for n in names}
    comp=prepare_phylogatr_compact_ttf_design(geos,train_species=train,eval_species=evaluation,bandwidth=500.,prior_strength=.25,segment_points=5,min_training_edges=5)
    tm={n:comp.template_edges[n].midpoint for n in train}; em={n:comp.template_edges[n].midpoint for n in evaluation}
    pools=prepare_target_source_pools(comp.prepared,tm,em,support_radius=500.,minimum_target_coverage=.5,minimum_source_species=5)
    kt,axes=trait_kernel(names,tr); kg,gfeat=geometry_kernel(names,geos); idx={n:i for i,n in enumerate(names)}
    pairs=[]
    for target in pools.eligible_eval_species:
        tc=np.mean(geos[target].coordinates,axis=0)
        for source in pools.source_pool[target]:
            sc=np.mean(geos[source].coordinates,axis=0)
            A=em[target]; B=tm[source]; dist2=np.sum((A[:,None,:]-B[None,:,:])**2,axis=2)
            ib=np.argmin(dist2,axis=1); ia=np.arange(len(A)); coverage=float(np.mean(np.sqrt(np.min(dist2,axis=1))<=500.))
            pairs.append({"target":target,"source":source,"trait_similarity":float(kt[idx[target],idx[source]]),"coverage":coverage,"centroid_distance":float(np.linalg.norm(tc-sc)),"edge_count_ratio":float(geos[source].n_edges/geos[target].n_edges),"locality_count_ratio":float(geos[source].n_localities/geos[target].n_localities),"target_edge_index":ia.tolist(),"source_edge_index":ib.tolist()})
    out={"schema":"ttf_lepidoptera_trait_gradient_design_v0.1","species_order":list(names),"coordinates":{n:geos[n].coordinates.tolist() for n in names},"train_species":list(train),"eval_species":list(evaluation),"eligible_eval_species":list(pools.eligible_eval_species),"trait_kernel":kt.tolist(),"geometry_kernel":kg.tolist(),"trait_axes":axes,"geometry_features":gfeat.tolist(),"pairs":pairs,"provenance":{"prior_identity_route_species_excluded":len(exclusion["species"]),"additional_operator_exposure_excluded":list(exclusion.get("additional_operator_exposure_exclusion",[]))},"outcome_firewall":{"sequence_identity_used":False,"genetic_distance_used":False,"empirical_transfer_used":False}}
    a.output.parent.mkdir(parents=True,exist_ok=True); a.output.write_text(json.dumps(out,separators=(",",":"))+"\n")
    print(json.dumps({"species":len(names),"train":len(train),"eval":len(evaluation),"eligible_eval":len(pools.eligible_eval_species),"pairs":len(pairs)}))
if __name__=="__main__":main()
