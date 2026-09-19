#!/usr/bin/env python3
from __future__ import annotations
import argparse,json
from pathlib import Path
import numpy as np
from ttf.lepidoptera_trait_gradient import residualize_within_target,equal_target_gradient
from ttf.lepidoptera_trait_gradient_simulate import make_worlds

def wilson(k,n,z=1.959963984540054):
    p=k/n; d=1+z*z/n; c=(p+z*z/(2*n))/d; h=z*np.sqrt(p*(1-p)/n+z*z/(4*n*n))/d
    return float(c-h),float(c+h)

def score_pair_contributions(world, pairs):
    # Frozen qualification surface: each row supplies target/source, four geometry
    # covariates, trait similarity, and an edge-alignment operator precomputed
    # without empirical response. The operator maps source and target synthetic
    # edge responses to one source-specific transfer contribution.
    y=[]; s=[]; g=[]; t=[]
    for p in pairs:
        a=np.asarray(world.edge_response[p["target"]]); b=np.asarray(world.edge_response[p["source"]])
        ia=np.asarray(p["target_edge_index"],int); ib=np.asarray(p["source_edge_index"],int)
        y.append(float(np.corrcoef(a[ia],b[ib])[0,1]) if len(ia)>=3 else 0.)
        s.append(float(p["trait_similarity"])); g.append([float(p[k]) for k in ("coverage","centroid_distance","edge_count_ratio","locality_count_ratio")]); t.append(p["target"])
    r=residualize_within_target(np.asarray(s),np.asarray(g),np.asarray(t))
    return equal_target_gradient(np.asarray(y),r,np.asarray(t))[0]

def main():
    ap=argparse.ArgumentParser(); ap.add_argument("--design-json",type=Path,required=True); ap.add_argument("--cell",choices=["private","trait_gradient_positive","geometry_confounded_trap"],required=True); ap.add_argument("--start",type=int,default=0); ap.add_argument("--count",type=int,required=True); ap.add_argument("--output",type=Path,required=True); a=ap.parse_args()
    d=json.loads(a.design_json.read_text()); from ttf.genetic_geometry import prepare_density_scaled_genetic_geometry
    geos={n:prepare_density_scaled_genetic_geometry(np.asarray(x,float),neighbor_fraction=.15) for n,x in d["coordinates"].items()}
    worlds=make_worlds(geos,d["trait_score"],d["geometry_score"],cell=a.cell,start=a.start,count=a.count)
    stats=np.asarray([score_pair_contributions(w,d["pairs"]) for w in worlds])
    out={"schema":"ttf_lepidoptera_trait_gradient_qualification_shard_v0.1","cell":a.cell,"start":a.start,"count":a.count,"statistics":stats.tolist(),"outcome_firewall":{"empirical_sequence_identity_opened":False,"empirical_transfer_statistic_computed":False}}
    a.output.parent.mkdir(parents=True,exist_ok=True); a.output.write_text(json.dumps(out,indent=2)+"\n")
if __name__=="__main__":main()
