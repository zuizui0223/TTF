#!/usr/bin/env python3
"""Synthetic proof that straight-line climate averages cannot determine corridors.

This is NOT a paleoclimate prediction, a host plant SDM, or a genetic test.
Four-neighbor steps and grid units are illustrative only; no real source,
parameterization or geographic domain is changed by this demonstration.
"""
from __future__ import annotations
import argparse
import heapq
import json
from pathlib import Path
import numpy as np


def illustrative_four_neighbor_least_cost(costs, start, end):
    """Illustrative Dijkstra on positive pixel costs; inf = impassable cell."""
    z=np.asarray(costs,dtype=float)
    if z.ndim!=2 or min(z.shape)<2 or np.isnan(z).any() or np.any(z<=0):
        raise ValueError("all traversable cells must have positive finite resistance")
    rows,cols=z.shape
    for rc in (start,end):
        if len(rc)!=2 or not (0<=rc[0]<rows and 0<=rc[1]<cols):
            raise ValueError("endpoint outside illustrative grid")
        if not np.isfinite(z[rc]):
            raise ValueError("endpoint is unreachable")
    dist=np.full(z.shape,np.inf)
    dist[start]=0
    pq=[(0.0,start[0],start[1])]
    while pq:
        d,r,c=heapq.heappop(pq)
        if d!=dist[r,c]:continue
        if (r,c)==tuple(end):return float(d)
        for dr,dc in ((0,-1),(0,1),(-1,0),(1,0)):
            rn,cn=r+dr,c+dc
            if not(0<=rn<rows and 0<=cn<cols) or not np.isfinite(z[rn,cn]):continue
            cand=d+0.5*(z[r,c]+z[rn,cn])
            if cand<dist[rn,cn]:
                dist[rn,cn]=cand
                heapq.heappush(pq,(cand,rn,cn))
    raise ValueError("unreachable endpoint under illustrative resistance mask")


def five_straight_line_samples(costs,start,end):
    """2D toy analogue of existing five interior edge fractions, NOT 2D GIS."""
    f=(0.1,0.3,0.5,0.7,0.9)
    c=np.asarray(costs,float)
    indexes=[(int(round((1-t)*start[0]+t*end[0])),
              int(round((1-t)*start[1]+t*end[1]))) for t in f]
    return float(np.mean([c[i] for i in indexes]))


def toy_counterexample():
    start,end=(4,0),(4,8)
    bad=np.full((9,9),20.0)
    bad[start]=bad[end]=1.0
    good=bad.copy()
    # A low-cost arc outside the straight segment exists ONLY in the first map.
    good[3,:]=1.0
    cost_good=five_straight_line_samples(good,start,end)
    cost_bad=five_straight_line_samples(bad,start,end)
    d_good=illustrative_four_neighbor_least_cost(good,start,end)
    d_bad=illustrative_four_neighbor_least_cost(bad,start,end)
    if cost_good!=cost_bad or not d_good<d_bad:
        raise RuntimeError("synthetic identifiability counterexample collapsed")
    return {
        "class":"SYNTHETIC_NONIDENTIFIABILITY_EXAMPLE",
        "five_point_straight_line_cost":{
            "with_alternative_corridor":cost_good,
            "without_alternative_corridor":cost_bad,
        },
        "least_cost_distance":{
            "with_alternative_corridor":d_good,
            "without_alternative_corridor":d_bad,
        },
        "units":"arbitrary artificial grid resistance per cardinal step",
        "not_an_empirical_LGM_result":True,
        "not_an_ancient_host_occurrence_model":True,
        "not_a_new_genetic_test":True,
        "illustrative_four_neighbor_only":True,
        "conclusion":"Equal five-site straight-line mean resistance does not identify connectivity when off-axis routes differ."
    }


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--output",type=Path,required=True)
    a=p.parse_args()
    result=toy_counterexample()
    a.output.parent.mkdir(parents=True,exist_ok=True)
    a.output.write_text(json.dumps(result,sort_keys=True,indent=2)+"\n",encoding="utf-8")
    print(json.dumps(result,sort_keys=True))

if __name__=="__main__":main()
