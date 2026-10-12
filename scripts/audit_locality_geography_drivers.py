#!/usr/bin/env python3
"""Response-exposed descriptive geographic sampling audit of locality deletions.

Joins prior (derived) site-deletion outcomes to untouched response-blind ECEF
locality/edge geometry. No nucleotide identities, gene distances or new tests.
The original 259-species confirmatory result stays NOT_EVALUABLE.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import statistics
from collections import defaultdict
from pathlib import Path

import numpy as np

SPECIES_STABILITY_SHA256 = "b447b360edc3d73ae62fe4c668f2cadf9c13f94a0292d8aee5e96828c6724f23"
SITE_DELETIONS_SHA256 = "ba6f0f5f1ec9b9da4e697a0a271f0156ac2ad80d385da0a793bba8e200f97686"
LOCATIONS_SHA256 = "037cd8fa1f059fb67c349a465540d3d5fac469b5d14a2a4d2658d74c036c0ae9"
EDGES_SHA256 = "ab10a876895cf00817e8ce555665ebb78a0f2ac64323d2e93c213e1857738ba9"

METRICS = ("nearest_neighbor_chord_km", "centroid_chord_km",
           "median_to_other_chord_km", "frozen_edge_degree")


def sha256_file(path:Path)->str:
    h=hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda:f.read(1<<20), b""):
            h.update(chunk)
    return h.hexdigest()


def rank_avg(values):
    """Average-tie rank / n within ONE species (not geographic distance units)."""
    v=np.asarray(values,dtype=float)
    if v.ndim!=1 or not len(v) or not np.isfinite(v).all():
        raise ValueError("nonfinite or empty frozen locality geography")
    order=np.argsort(v,kind="mergesort")
    ranks=np.empty(len(v),float)
    i=0
    while i<len(v):
        j=i+1
        while j<len(v) and v[order[i]]==v[order[j]]:
            j+=1
        ranks[order[i:j]]=(i+1+j)/2/len(v)
        i=j
    return ranks


def species_geographic_ranked_features(coords, pairs):
    xyz=np.asarray(coords,dtype=float)
    e=np.asarray(pairs,dtype=int)
    n=len(xyz)
    if (xyz.shape!=(n,3) or e.ndim!=2 or e.shape[1]!=2
            or n<2 or not np.isfinite(xyz).all()):
        raise ValueError("invalid frozen ECEF geometry")
    if np.any(e<0) or np.any(e>=n):
        raise ValueError("edge points to nonexistent frozen locality")
    if np.any(e[:,0]==e[:,1]):
        raise ValueError("frozen graph self loop")
    degree=np.bincount(e.ravel(),minlength=n).astype(float)
    dist=np.linalg.norm(xyz[:,None,:]-xyz[None,:,:],axis=2)
    np.fill_diagonal(dist,np.inf)
    near=dist.min(axis=1)
    # Use all other distances; infinity is the excluded diagonal only.
    distant=np.where(np.isinf(dist),np.nan,dist)
    median=np.nanmedian(distant,axis=1)
    to_centroid=np.linalg.norm(xyz-xyz.mean(axis=0),axis=1)
    return {
        "nearest_neighbor_chord_km":rank_avg(near),
        "centroid_chord_km":rank_avg(to_centroid),
        "median_to_other_chord_km":rank_avg(median),
        "frozen_edge_degree":rank_avg(degree)
    }


def analyze_geographic_influence(species_rows, deletion_rows, localities, edges):
    sample={}
    for s in species_rows:
        key=s["species"]
        if key in sample: raise ValueError("duplicate species in frozen response-exposed coefficients")
        sample[key]=s
    if not sample: raise ValueError("empty frozen confirmatory set")

    # Full 642-source geometry may include non-survivors: select only exact
    # surviving 259 species without ever reading their nucleotide identities.
    loc=defaultdict(dict)
    for l in localities:
        sp=l["species"]
        if sp not in sample: continue
        idx=int(l["locality_index"])
        if idx in loc[sp]: raise ValueError("duplicate frozen locality")
        loc[sp][idx]=[float(l[k]) for k in ("x_km","y_km","z_km")]
    edge=defaultdict(dict)
    for row in edges:
        sp=row["species"]
        if sp not in sample: continue
        idx=int(row["edge_index"])
        if idx in edge[sp]: raise ValueError("duplicate frozen edge")
        edge[sp][idx]=(int(row["node_left"]),int(row["node_right"]))
    by_delete=defaultdict(dict)
    for d in deletion_rows:
        sp=d["species"]
        if sp not in sample:
            raise ValueError("deletion species absent from original survivor coefficients")
        idx=int(d["locality_index"])
        if idx in by_delete[sp]: raise ValueError("duplicate deletion locality")
        by_delete[sp][idx]=d

    metric_flip=defaultdict(list)
    metric_other=defaultdict(list)
    paired=defaultdict(list)
    zero=0
    flips=0
    flipped_species=0
    defined=0
    undefined=0
    total_edges=0

    for sp,row in sample.items():
        n=int(row["localities"])
        e_n=int(row["edges"])
        if (set(loc.get(sp,{}))!=set(range(n))
                or set(by_delete.get(sp,{}))!=set(range(n))
                or set(edge.get(sp,{}))!=set(range(e_n))):
            raise ValueError(f"frozen locality or edge count/index mismatch: {sp}")
        coord=[loc[sp][i] for i in range(n)]
        e=[edge[sp][i] for i in range(e_n)]
        ranked=species_geographic_ranked_features(coord,e)
        base=float(row["unscaled_slope"])
        if not math.isfinite(base): raise ValueError("nonfinite original slope")
        if base==0: zero+=1
        else: total_edges+=e_n
        yes=[];no=[]
        for j in range(n):
            d=by_delete[sp][j]
            if d["status"]=="undefined":
                undefined+=1
                continue
            if d["status"]!="defined": raise ValueError("illegal delete status")
            v=float(d["slope"])
            if not math.isfinite(v): raise ValueError("nonfinite deletion coefficient")
            if base==0:
                if abs(v)>1e-13: raise ValueError("formerly invariant slope changed unexpectedly")
                continue
            defined+=1
            if base*v<0:
                yes.append(j);flips+=1
            else: no.append(j)
        if base!=0 and yes:
            flipped_species+=1
        if int(row["locality_deletions_reversing_sign"])!=len(yes):
            raise ValueError("frozen strict sign flips do not match source")
        for metric,rank in ranked.items():
            metric_flip[metric].extend(rank[yes].tolist())
            metric_other[metric].extend(rank[no].tolist())
            if yes and no:
                paired[metric].append(float(np.mean(rank[yes])-np.mean(rank[no])))

    summary={}
    for metric in METRICS:
        fs=metric_flip[metric];ns=metric_other[metric]
        paired_diffs=paired[metric]
        summary[metric]={
            "flip_site_n":len(fs),"other_site_n":len(ns),
            "flipping_site_mean_rank":float(statistics.mean(fs)) if fs else None,
            "other_site_mean_rank":float(statistics.mean(ns)) if ns else None,
            "paired_species":len(paired_diffs),
            "paired_median_rank_difference":float(statistics.median(paired_diffs)) if paired_diffs else None,
            "species_with_higher_flip_site_rank":sum(x>0 for x in paired_diffs),
            "species_with_lower_flip_site_rank":sum(x<0 for x in paired_diffs),
            "species_with_tied_flip_site_rank":sum(x==0 for x in paired_diffs),
        }
    return {
        "schema":"ttf_hhm259_post_outcome_geographic_locality_influence_v0.1",
        "status":"EXPLORATORY_SAMPLING_GEOGRAPHY_NOT_HISTORICAL_BARRIER_TEST",
        "original_confirmatory_status":"NOT_EVALUABLE_EMPIRICAL_RESPONSE",
        "original_species":len(sample),
        "original_zero_species":zero,
        "original_sign_flipping_species":flipped_species,
        "site_deletions_reversing_sign":flips,
        "defined_deletions_nonzero":defined,
        "undefined_deletions":undefined,
        "within_species_geographic_comparisons":summary,
        "interpretation":{
            "nearest_chord_rank_proxy_for_local_sampling_isolation":True,
            "different_species_geographic_range_normalized_by_rank":True,
            "site_results_are_clustered_within_species_not_independent_replicates":True,
            "no_p_values_or_new_significance_claims":True,
            "does_not_demonstrate_paleogeographic_barriers":True,
            "all_frozen_259_species_and_edges_preserved":True,
            "no_genetic_response_reopening":True
        }
    }


def main():
    a=argparse.ArgumentParser()
    for name in ("species_stability","site_deletions","frozen_localities","frozen_edges"):
        a.add_argument("--"+name.replace("_","-"),type=Path,required=True)
    a.add_argument("--output",type=Path,required=True)
    args=a.parse_args()
    pins={
        "species_stability":SPECIES_STABILITY_SHA256,
        "site_deletions":SITE_DELETIONS_SHA256,
        "frozen_localities":LOCATIONS_SHA256,
        "frozen_edges":EDGES_SHA256,
    }
    for key,sha in pins.items():
        if sha256_file(getattr(args,key))!=sha:
            raise RuntimeError(f"frozen input hash mismatch: {key}")
    def read_csv(path):
        with path.open(newline="",encoding="utf-8") as f:
            return list(csv.DictReader(f))
    summary=analyze_geographic_influence(
        read_csv(args.species_stability),read_csv(args.site_deletions),
        read_csv(args.frozen_localities),read_csv(args.frozen_edges)
    )
    if (summary["original_species"],summary["original_zero_species"],
            summary["original_sign_flipping_species"],
            summary["site_deletions_reversing_sign"],
            summary["defined_deletions_nonzero"],
            summary["undefined_deletions"])!=(259,4,129,420,6308,1):
        raise RuntimeError("source-bound exposed geography audit drift")
    summary["input_sha256"]=pins
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(summary,sort_keys=True,indent=2)+"\n")
    print(json.dumps({
        "status":summary["status"],
        "n_flip_species":summary["original_sign_flipping_species"],
        "nearest_neighbor":summary["within_species_geographic_comparisons"]["nearest_neighbor_chord_km"]
    },sort_keys=True))

if __name__ == "__main__":
    main()
