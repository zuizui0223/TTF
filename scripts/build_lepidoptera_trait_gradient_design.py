#!/usr/bin/env python3
from __future__ import annotations
import argparse,csv,hashlib,json,zipfile,tempfile
from pathlib import Path
import numpy as np
from ttf.phylogatr_confirmatory import scan_phylogatr_phase1,choose_one_panel_per_species
from ttf.genetic_conditional_geometry_match import prepare_conditional_geometry_matched_order
from ttf.lepidoptera_trait_gradient import composite_trait_similarity

ALIASES=("COI","CO1","COX1","COXI","CYTOCHROME C OXIDASE SUBUNIT I","CYTOCHROME C OXIDASE SUBUNIT 1")
EXCLUDE={"Acrolophus arcanella"}

def traits(path):
    out={}
    with path.open(newline="",encoding="utf-8-sig") as f:
        for x in csv.DictReader(f):
            n=x["Species"].strip()
            if n and n not in out: out[n]=x
    return out

def wing(x):
    vals=[x.get(k,"") for k in ("WS_L","WS_U","FW_L","FW_U")]
    z=[float(v) for v in vals if str(v).strip() not in {"","NA"}]
    return float(np.mean(z))
def habitat(x): return tuple(x[k] for k in ("CanopyAffinity","EdgeAffinity","MoistureAffinity","DisturbanceAffinity"))
def complete(x):
    try: wing(x); float(x["NumberOfHostplantFamilies"])
    except: return False
    return str(x["Voltinism"]).strip() not in {"","NA"} and all(str(v).strip() not in {"","NA"} for v in habitat(x))

def main():
    ap=argparse.ArgumentParser(); ap.add_argument("--archive",type=Path,required=True); ap.add_argument("--leptraits",type=Path,required=True); ap.add_argument("--output",type=Path,required=True); a=ap.parse_args()
    with tempfile.TemporaryDirectory() as td:
        with zipfile.ZipFile(a.archive) as z:z.extractall(td)
        roots=[p.parent for p in Path(td).rglob("genes.txt") if (p.parent/"cite.txt").is_file()]
        if len(roots)!=1: raise RuntimeError("archive root drift")
        scan=scan_phylogatr_phase1(roots[0],aliases=ALIASES,excluded_species=EXCLUDE,min_localities=12,min_endpoint_training_edges=5,neighbor_fraction=.15)
        panels=choose_one_panel_per_species(scan.candidates)
    tr=traits(a.leptraits); chosen={p.species:p for p in panels if p.order=="Lepidoptera" and p.species in tr and complete(tr[p.species])}
    names=tuple(sorted(chosen,key=lambda n:(hashlib.sha256(f"lepidoptera-trait-transfer-v0.1|major_complete|{n}".encode()).hexdigest(),n)))
    cut=len(names)//2; train=names[:cut]; evaluation=names[cut:]
    coords={n:chosen[n].geometry.coordinates.tolist() for n in names}
    wingv={n:wing(tr[n]) for n in names}; vol={n:tr[n]["Voltinism"] for n in names}; host={n:float(tr[n]["NumberOfHostplantFamilies"]) for n in names}; hab={n:habitat(tr[n]) for n in names}
    # geometry support is encoded pairwise without response.
    from ttf.conditional_transfer import prepare_target_source_pools
    from ttf.phylogatr_compact_execution import prepare_phylogatr_compact_ttf_design
    geos={n:chosen[n].geometry for n in names}; comp=prepare_phylogatr_compact_ttf_design(geos,train_species=train,eval_species=evaluation,bandwidth=500.,prior_strength=.25,segment_points=5,min_training_edges=5)
    tm={n:comp.template_edges[n].midpoint for n in train}; em={n:comp.template_edges[n].midpoint for n in evaluation}
    pools=prepare_target_source_pools(comp.prepared,tm,em,support_radius=500.,minimum_target_coverage=.5,minimum_source_species=5)
    pairs=[]
    for target in pools.eligible_eval_species:
        tc=np.mean(geos[target].coordinates,axis=0)
        for source in pools.source_pool[target]:
            sc=np.mean(geos[source].coordinates,axis=0)
            # nearest-edge alignment, response blind
            A=em[target]; B=tm[source]; dist=np.sum((A[:,None,:]-B[None,:,:])**2,axis=2); ib=np.argmin(dist,axis=1); ia=np.arange(len(A))
            coverage=float(np.mean(np.sqrt(np.min(dist,axis=1))<=500.))
            pairs.append({"target":target,"source":source,"trait_similarity":composite_trait_similarity(target,source,wing_size=wingv,voltinism=vol,host_breadth=host,habitat=hab),"coverage":coverage,"centroid_distance":float(np.linalg.norm(tc-sc)),"edge_count_ratio":float(geos[source].n_edges/geos[target].n_edges),"locality_count_ratio":float(geos[source].n_localities/geos[target].n_localities),"target_edge_index":ia.tolist(),"source_edge_index":ib.tolist()})
    # species-level simulator scores are fixed response-blind summaries.
    ts={n:float(np.mean([p["trait_similarity"] for p in pairs if p["source"]==n] or [0.])) for n in names}
    gs={n:float(np.mean([np.exp(-p["centroid_distance"]/500.) for p in pairs if p["source"]==n] or [0.])) for n in names}
    out={"schema":"ttf_lepidoptera_trait_gradient_design_v0.1","coordinates":coords,"train_species":list(train),"eval_species":list(evaluation),"eligible_eval_species":list(pools.eligible_eval_species),"trait_score":ts,"geometry_score":gs,"pairs":pairs,"outcome_firewall":{"sequence_identity_used":False,"genetic_distance_used":False,"empirical_transfer_used":False}}
    a.output.parent.mkdir(parents=True,exist_ok=True); a.output.write_text(json.dumps(out,separators=(",",":"))+"\n")
    print(json.dumps({"species":len(names),"train":len(train),"eval":len(evaluation),"eligible_eval":len(pools.eligible_eval_species),"pairs":len(pairs)}))
if __name__=="__main__":main()
