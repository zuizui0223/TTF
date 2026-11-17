"""Design-only precision for restoring Ellis 2007 within-basin cross-quartz arms.

Exact Fraction arithmetic; iid homoscedastic unit-error working model only.
No fitness/genetics/GPS/LGM values opened.
"""
from __future__ import annotations
import argparse,hashlib,itertools,json
from fractions import Fraction
from pathlib import Path
from ellis2007_design_rank_v01 import RECEIPT,SHA,analyze,category,CLASSES

def inverse(a):
    n=len(a)
    m=[[Fraction(v) for v in row]+[Fraction(int(i==j)) for j in range(n)] for i,row in enumerate(a)]
    for j in range(n):
        p=next((i for i in range(j,n) if m[i][j]),None)
        if p is None:raise ValueError("singular design")
        m[j],m[p]=m[p],m[j]
        z=m[j][j];m[j]=[x/z for x in m[j]]
        for i in range(n):
            if i!=j:
                z=m[i][j];m[i]=[x-z*y for x,y in zip(m[i],m[j])]
    return [row[n:] for row in m]

def row(s,g):
    return [1,s[1],g[1],int(s[0]==1),int(s[0]==2),int(g[0]==1),int(g[0]==2)]

def variance(edges):
    x=[row(*edge) for edge in edges]
    gram=[[sum(a[i]*a[j] for a in x) for j in range(7)] for i in range(7)]
    inv=inverse(gram)
    c=[0,1,-1,0,0,0,0]
    return sum(Fraction(c[i])*inv[i][j]*c[j] for i in range(7) for j in range(7))

def evaluate(data):
    analyze(data)
    nodes=list(itertools.product(range(3),range(2)))
    base=[(s,g) for s in nodes for g in nodes if category(s,g) in (CLASSES[0],CLASSES[2])]
    missing=[(s,g) for s in nodes for g in nodes if category(s,g)==CLASSES[1]]
    if len(base)!=18 or len(missing)!=6:raise AssertionError("frozen topology changed")
    grouped={}
    for k in range(1,7):
        values=sorted(variance(base+list(extra)) for extra in itertools.combinations(missing,k))
        grouped[str(k)]={
          "subsets":len(values),
          "min_variance_factor":float(values[0]),
          "median_variance_factor":float(values[len(values)//2]) if len(values)%2 else float((values[len(values)//2-1]+values[len(values)//2])/2),
          "max_variance_factor":float(values[-1])
        }
    if grouped["1"]["subsets"]!=6 or grouped["6"]["subsets"]!=1:raise AssertionError("subset count drift")
    if variance(base+missing)!=Fraction(2,3):raise AssertionError("full design variance drift")
    return {
      "schema":"ttf_Ellis2007_source_only_quartz_precision_v01",
      "source_receipt_sha256":SHA,
      "status":"DESIGN_ONLY_WORKING_MODEL_VARIANCE_NOT_EMPIRICAL_PRECISION",
      "contrast":"source_quartz_minus_garden_quartz",
      "working_assumption":"independent homoscedastic errors with common variance sigma_squared; variance_factor multiplies sigma_squared",
      "missing_arm_subsets_exhausted":63,
      "baseline_zero_restored_arms":"NONESTIMABLE",
      "by_restored_arm_count":grouped,
      "one_arm_vs_all_six_variance_factor_ratio":float(variance(base+[missing[0]])/variance(base+missing)),
      "original_fitness_genetics_GPS_LGM_outcomes_opened":False,
      "independent_replicate_or_causal_precision_claim":False,
      "frozen_endpoints_unchanged":True
    }

def main():
    p=argparse.ArgumentParser();p.add_argument("--output",required=True,type=Path);a=p.parse_args()
    raw=RECEIPT.read_bytes()
    if hashlib.sha256(raw).hexdigest()!=SHA:raise ValueError("receipt SHA drift")
    r=evaluate(json.loads(raw));a.output.parent.mkdir(parents=True,exist_ok=True)
    a.output.write_text(json.dumps(r,sort_keys=True,indent=2)+"\n",encoding="utf-8")
    print(json.dumps(r,sort_keys=True))
if __name__=="__main__":main()
