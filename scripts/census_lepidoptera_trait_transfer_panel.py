#!/usr/bin/env python3
from __future__ import annotations
import argparse,csv,hashlib,json
from pathlib import Path
import numpy as np
from ttf.conditional_transfer import prepare_target_source_pools
from ttf.genetic_geometry import prepare_density_scaled_genetic_geometry
from ttf.phylogatr_compact_execution import prepare_phylogatr_compact_ttf_design

SCHEMA="ttf_lepidoptera_trait_transfer_census_v0.1"
MISS={"","NA","N/A","NULL","NONE","NAN"}
HAB=("CanopyAffinity","EdgeAffinity","MoistureAffinity","DisturbanceAffinity")
WING=("WS_L","WS_U","FW_L","FW_U","WS_L_Fem","WS_U_Fem","WS_L_Mal","WS_U_Mal","FW_L_Fem","FW_U_Fem","FW_L_Mal","FW_U_Mal")

def present(x): return str(x).strip().upper() not in MISS
def sha(path): return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def load_rows(path):
    out={}
    with Path(path).open(newline="",encoding="utf-8") as f:
        r=csv.DictReader(f); req={"species","order","locality_index","x_km","y_km","z_km","graph_k","edge_count","min_endpoint_disjoint_training_edges"}
        if not req<=set(r.fieldnames or ()): raise RuntimeError("geometry CSV schema mismatch")
        for x in r: out.setdefault(x["species"].strip(),[]).append(x)
    return out

def eligible_leps(rows,min_localities=12,min_edges=5):
    ans=[]
    for sp,rs in rows.items():
        if {x["order"].strip() for x in rs}!={"Lepidoptera"}: continue
        if len({int(x["locality_index"]) for x in rs})<min_localities: continue
        if min(int(x["min_endpoint_disjoint_training_edges"]) for x in rs)<min_edges: continue
        ans.append(sp)
    return tuple(sorted(ans))

def load_traits(path):
    out={}
    with Path(path).open(newline="",encoding="utf-8-sig") as f:
        r=csv.DictReader(f)
        if "Species" not in set(r.fieldnames or ()): raise RuntimeError("LepTraits CSV missing Species")
        for x in r:
            n=x["Species"].strip()
            if n and n not in out: out[n]=x
    return out

def load_map(path):
    if path is None:return {}
    out={}
    with Path(path).open(newline="",encoding="utf-8") as f:
        r=csv.DictReader(f)
        if not {"source_species","canonical_species"}<=set(r.fieldnames or ()): raise RuntimeError("name map schema mismatch")
        for x in r:
            if x["source_species"].strip() and x["canonical_species"].strip(): out[x["source_species"].strip()]=x["canonical_species"].strip()
    return out

def flags(x):
    wing=any(present(x.get(k,"")) for k in WING); vol=present(x.get("Voltinism","")); host=present(x.get("NumberOfHostplantFamilies",""))
    h=[present(x.get(k,"")) for k in HAB]
    return {"wing_size":wing,"voltinism":vol,"host_breadth":host,"habitat_any":any(h),"habitat_complete4":all(h),"major_complete":wing and vol and host and all(h)}

def split(names,tag):
    z=sorted(names,key=lambda n:(hashlib.sha256(f"{tag}|{n}".encode()).hexdigest(),n)); c=len(z)//2
    return tuple(z[:c]),tuple(z[c:])

def geos(rows,names):
    out={}
    for n in names:
        rs=sorted(rows[n],key=lambda x:int(x["locality_index"]))
        xyz=np.asarray([[float(x["x_km"]),float(x["y_km"]),float(x["z_km"])] for x in rs])
        g=prepare_density_scaled_genetic_geometry(xyz,neighbor_fraction=.15)
        if {g.graph_k}!={int(x["graph_k"]) for x in rs} or {g.n_edges}!={int(x["edge_count"]) for x in rs}: raise RuntimeError(f"geometry drift: {n}")
        out[n]=g
    return out

def support(rows,names,tag,radius=500.,coverage=.5,min_sources=5):
    if len(names)<4:return {"train_species":0,"eval_species":0,"supported_eval_species":0,"source_count_quantiles":None}
    tr,ev=split(tuple(names),tag); comp=prepare_phylogatr_compact_ttf_design(geos(rows,tuple(names)),train_species=tr,eval_species=ev,bandwidth=radius,prior_strength=.25,segment_points=5,min_training_edges=5)
    tm={n:comp.template_edges[n].midpoint for n in comp.train_species}; em={n:comp.template_edges[n].midpoint for n in comp.eval_species}
    p=prepare_target_source_pools(comp.prepared,tm,em,support_radius=radius,minimum_target_coverage=coverage,minimum_source_species=min_sources)
    a=np.asarray([len(p.source_pool[n]) for n in comp.eval_species],float)
    return {"train_species":len(tr),"eval_species":len(ev),"supported_eval_species":len(p.eligible_eval_species),"source_count_quantiles":dict(zip(("min","q25","median","q75","max"),map(float,np.quantile(a,[0,.25,.5,.75,1]))))}

def main():
    q=argparse.ArgumentParser(); q.add_argument("--geometry-csv",type=Path,required=True); q.add_argument("--leptraits-csv",type=Path,required=True); q.add_argument("--name-map-csv",type=Path); q.add_argument("--output",type=Path,required=True); q.add_argument("--split-tag",default="lepidoptera-trait-transfer-v0.1"); a=q.parse_args()
    rows=load_rows(a.geometry_csv); leps=eligible_leps(rows); traits=load_traits(a.leptraits_csv); nm=load_map(a.name_map_csv)
    matched={n:flags(traits[nm.get(n,n)]) for n in leps if nm.get(n,n) in traits}
    stages={"metadata_eligible_lepidoptera":leps,"leptraits_name_matched":tuple(sorted(matched))}
    for k in ("wing_size","voltinism","host_breadth","habitat_any","habitat_complete4","major_complete"): stages[k]=tuple(sorted(n for n,v in matched.items() if v[k]))
    sup={k:support(rows,stages[k],f"{a.split_tag}|{k}") for k in ("leptraits_name_matched","wing_size","voltinism","host_breadth","habitat_complete4","major_complete")}
    n=sup["major_complete"]["supported_eval_species"]; decision="GREEN_GE_90" if n>=90 else "AMBER_60_TO_89" if n>=60 else "RED_LE_40" if n<=40 else "INTERMEDIATE_41_TO_59"
    out={"schema":SCHEMA,"status":"RESPONSE_BLIND_CENSUS_ONLY","inputs":{"geometry_csv_sha256":sha(a.geometry_csv),"leptraits_csv_sha256":sha(a.leptraits_csv),"name_map_csv_sha256":None if a.name_map_csv is None else sha(a.name_map_csv)},"funnel_counts":{k:len(v) for k,v in stages.items()},"geographic_support":sup,"major_complete_decision":decision,"outcome_firewall":{"sequence_identity_opened":False,"pairwise_genetic_distances_opened":False,"synthetic_response_worlds_opened":False,"empirical_transfer_statistic_computed":False}}
    a.output.parent.mkdir(parents=True,exist_ok=True); a.output.write_text(json.dumps(out,indent=2,sort_keys=True)+"\n"); print(json.dumps({"funnel_counts":out["funnel_counts"],"decision":decision},indent=2))
if __name__=="__main__": main()
