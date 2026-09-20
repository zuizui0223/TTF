#!/usr/bin/env python3
from __future__ import annotations
import argparse,hashlib,json
from pathlib import Path
import numpy as np
from ttf.core import average_ranks
from ttf.genetic_geometry import prepare_density_scaled_genetic_geometry
from ttf.lepidoptera_trait_gradient import residualize_within_target,equal_target_gradient
from ttf.lepidoptera_trait_gradient_simulate import prepare_trait_gradient_simulator,simulate_prepared_trait_gradient_world,frozen_seed

def wilson(k,n,z=1.959963984540054):
    p=k/n; d=1+z*z/n; c=(p+z*z/(2*n))/d; h=z*np.sqrt(p*(1-p)/n+z*z/(4*n*n))/d
    return float(c-h),float(c+h)

def bootstrap_seed(cell,replicate):
    x=f"20260920|lepidoptera-trait-gradient-v01-bootstrap|{cell}|{replicate}".encode()
    return int.from_bytes(hashlib.sha256(x).digest()[:8],"big")

def centered_two_sided_target_bootstrap(slopes,*,seed,n_bootstrap=1999):
    x=np.asarray(slopes,float); x=x[np.isfinite(x)]
    if len(x)<6: raise ValueError("at least six target slopes required")
    mean=float(np.mean(x)); sd=float(np.std(x,ddof=1)); se=sd/np.sqrt(len(x))
    obs=0.0 if se<=np.finfo(float).tiny and mean==0 else (np.inf*np.sign(mean) if se<=np.finfo(float).tiny else mean/se)
    centered=x-mean; rng=np.random.default_rng(int(seed)); ind=rng.integers(0,len(x),size=(int(n_bootstrap),len(x)))
    draws=centered[ind]; means=np.mean(draws,axis=1); sds=np.std(draws,axis=1,ddof=1); ses=sds/np.sqrt(len(x))
    t=np.zeros(int(n_bootstrap)); ok=ses>np.finfo(float).tiny; t[ok]=means[ok]/ses[ok]
    t[~ok & (means>0)]=np.inf; t[~ok & (means<0)]=-np.inf
    p=float((1+np.count_nonzero(np.abs(t)>=abs(obs)))/(int(n_bootstrap)+1))
    return mean,p

def prepare_pair_surface(pairs):
    similarity=np.asarray([float(p["trait_similarity"]) for p in pairs])
    geometry=np.asarray([[float(p[k]) for k in ("coverage","centroid_distance","edge_count_ratio","locality_count_ratio")] for p in pairs])
    target=np.asarray([str(p["target"]) for p in pairs])
    residual=residualize_within_target(similarity,geometry,target)
    pair_target=tuple(str(p["target"]) for p in pairs)
    pair_source=tuple(str(p["source"]) for p in pairs)
    source_index=tuple(np.asarray(p["source_edge_index"],dtype=np.int64) for p in pairs)
    return target,residual,pair_target,pair_source,source_index

def score_pair_contributions(world,pairs,target,residual,pair_target,pair_source,source_index):
    target_rank={}
    for name in set(pair_target):
        rank=np.asarray(average_ranks(np.asarray(world.edge_response[name])),dtype=float)
        centered=rank-float(rank.mean())
        target_rank[name]=(centered,float(np.dot(centered,centered)))
    y=np.empty(len(pairs),float)
    for i in range(len(pairs)):
        dx,xx=target_rank[pair_target[i]]
        aligned=np.asarray(world.edge_response[pair_source[i]])[source_index[i]]
        if len(aligned)!=len(dx): raise RuntimeError("nearest-edge alignment length drift")
        if len(aligned)<3:
            y[i]=0.0
            continue
        rank=np.asarray(average_ranks(aligned),dtype=float)
        dy=rank-float(rank.mean())
        den=float(np.sqrt(xx*np.dot(dy,dy)))
        y[i]=0.0 if den<=np.finfo(float).eps else float(np.dot(dx,dy)/den)
    stat,slopes=equal_target_gradient(y,residual,target)
    return float(stat),np.asarray(list(slopes.values()),float)

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--design-json",type=Path,required=True)
    ap.add_argument("--cell",choices=["private","trait_gradient_positive","geometry_confounded_trap"],required=True)
    ap.add_argument("--start",type=int,default=0); ap.add_argument("--count",type=int,required=True)
    ap.add_argument("--output",type=Path,required=True); a=ap.parse_args()
    d=json.loads(a.design_json.read_text())
    if d.get("schema")!="ttf_lepidoptera_trait_gradient_design_v0.1": raise RuntimeError("design schema drift")
    if any(d["outcome_firewall"].values()): raise RuntimeError("design firewall open")
    geos={n:prepare_density_scaled_genetic_geometry(np.asarray(x,float),neighbor_fraction=.15) for n,x in d["coordinates"].items()}
    species_order=tuple(map(str,d["species_order"])); tk=np.asarray(d["trait_kernel"],float); gk=np.asarray(d["geometry_kernel"],float)
    target,residual,pair_target,pair_source,source_index=prepare_pair_surface(d["pairs"])
    simulator=prepare_trait_gradient_simulator(geos,species_order,tk,gk,shared_fraction=0.85)
    stats=[]; pvals=[]; target_counts=[]
    for offset in range(a.count):
        replicate=a.start+offset
        w=simulate_prepared_trait_gradient_world(
            simulator,cell=a.cell,seed=frozen_seed(20260920,a.cell,replicate),
            private_amplitude=0.35,noise_sd=0.10,transition_width=0.20,latent_fields=6,
        )
        stat,slopes=score_pair_contributions(
            w,d["pairs"],target,residual,pair_target,pair_source,source_index
        )
        mean,p=centered_two_sided_target_bootstrap(
            slopes,seed=bootstrap_seed(a.cell,replicate),n_bootstrap=1999
        )
        if not np.isclose(stat,mean,atol=1e-12,rtol=0): raise RuntimeError("target slope mean drift")
        stats.append(stat); pvals.append(p); target_counts.append(len(slopes))
    stats=np.asarray(stats); pvals=np.asarray(pvals)
    if a.cell=="trait_gradient_positive":
        rejected=(pvals<=0.05)&(stats>0)
    else:
        rejected=pvals<=0.05
    k=int(np.count_nonzero(rejected)); lo,hi=wilson(k,len(rejected))
    gate=None
    if a.start==0 and a.count==500:
        gate=(lo>=0.80) if a.cell=="trait_gradient_positive" else (hi<=0.10)
    out={
        "schema":"ttf_lepidoptera_trait_gradient_qualification_shard_v0.1",
        "cell":a.cell,"start":a.start,"count":a.count,
        "statistics":stats.tolist(),"p_values":pvals.tolist(),
        "target_count_min":int(min(target_counts)),"target_count_max":int(max(target_counts)),
        "rejections":k,"rejection_rate":float(k/len(rejected)),
        "wilson95_lower":lo,"wilson95_upper":hi,"gate_pass":gate,
        "bootstrap":{"resamples":1999,"two_sided":True,"sampling_unit":"evaluation target"},
        "outcome_firewall":{"empirical_sequence_identity_opened":False,"empirical_transfer_statistic_computed":False},
    }
    a.output.parent.mkdir(parents=True,exist_ok=True)
    a.output.write_text(json.dumps(out,indent=2)+"\n")
    print(json.dumps({k:out[k] for k in ("cell","rejections","rejection_rate","wilson95_lower","wilson95_upper","gate_pass")},indent=2))
if __name__=="__main__": main()
