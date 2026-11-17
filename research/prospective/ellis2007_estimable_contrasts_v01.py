"""Exact-rational estimable contrasts for Ellis 2007 transplant design.

Only frozen source-only topology is accessed. No fitness or historical response.
"""
from __future__ import annotations
import argparse,hashlib,json
from fractions import Fraction
from itertools import product
from pathlib import Path
from ellis2007_design_rank_v01 import RECEIPT,SHA,analyze,rank,category,CLASSES

HERE=Path(__file__).resolve().parent
COLUMN_NAMES=["intercept","source_quartz","garden_quartz","source_basin_2","source_basin_3","garden_basin_2","garden_basin_3"]

def design_row(s,g):
    return [1,s[1],g[1],int(s[0]==1),int(s[0]==2),int(g[0]==1),int(g[0]==2)]

def null_basis(a):
    m=[[Fraction(x) for x in row] for row in a]
    if not m:return []
    n=len(m[0]);r=0;pivots=[]
    for j in range(n):
        k=next((i for i in range(r,len(m)) if m[i][j]),None)
        if k is None:continue
        m[r],m[k]=m[k],m[r]
        z=m[r][j];m[r]=[v/z for v in m[r]]
        for i in range(len(m)):
            if i!=r and m[i][j]:
                z=m[i][j];m[i]=[v-z*w for v,w in zip(m[i],m[r])]
        pivots.append(j);r+=1
        if r==len(m):break
    basis=[]
    for free in range(n):
        if free in pivots:continue
        v=[Fraction(0) for _ in range(n)];v[free]=1
        for i,j in enumerate(pivots):v[j]=-m[i][free]
        basis.append(v)
    return basis

def estimable(c,basis):
    return all(sum(Fraction(x)*y for x,y in zip(c,b))==0 for b in basis)

def evaluate(receipt):
    analyze(receipt)
    nodes=list(product(range(3),range(2)))
    primary=[design_row(s,g) for s in nodes for g in nodes if category(s,g) in CLASSES[:3]]
    secondary=[design_row(s,g) for s in nodes for g in nodes if category(s,g) in (CLASSES[0],CLASSES[2])]
    pb=null_basis(primary);sb=null_basis(secondary)
    if len(pb)!=0 or len(sb)!=1 or sb[0]!=[0,-1,1,0,0,0,0]:
        raise AssertionError("unexpected null space")
    contrasts={
      "source_quartz_only":[0,1,0,0,0,0,0],
      "garden_quartz_only":[0,0,1,0,0,0,0],
      "source_plus_garden_quartz":[0,1,1,0,0,0,0],
      "source_minus_garden_quartz":[0,1,-1,0,0,0,0],
      "source_basin_2_vs_1":[0,0,0,1,0,0,0],
      "garden_basin_2_vs_1":[0,0,0,0,0,1,0],
    }
    return {
      "schema":"ttf_Ellis2007_source_only_estimable_contrasts_v01",
      "source_receipt_sha256":SHA,
      "status":"ALGEBRAIC_CONTRAST_IDENTIFIABILITY_NOT_CAUSAL_OR_FITNESS_INFERENCE",
      "columns":COLUMN_NAMES,
      "primary_rank":rank(primary),"secondary_rank":rank(secondary),
      "primary_nullity":len(pb),"secondary_nullity":len(sb),
      "secondary_null_vector":[int(x) for x in sb[0]],
      "contrasts":{k:{"primary_estimable":estimable(c,pb),"secondary_estimable":estimable(c,sb)}
                   for k,c in contrasts.items()},
      "all_original_outcomes_unopened":True,
      "no_new_climate_or_LGM_effect":True,
      "frozen_endpoints_unchanged":True,
      "interpretation":"Only the sum of separate source and garden quartz coefficients is estimable in the secondary additive design. Not an identified adaptation effect."
    }

def main():
    p=argparse.ArgumentParser();p.add_argument("--output",required=True,type=Path)
    args=p.parse_args()
    raw=RECEIPT.read_bytes()
    if hashlib.sha256(raw).hexdigest()!=SHA:raise ValueError("receipt SHA drift")
    r=evaluate(json.loads(raw))
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(r,sort_keys=True,indent=2)+"\n",encoding="utf-8")
    print(json.dumps(r,sort_keys=True))
if __name__=="__main__":main()
