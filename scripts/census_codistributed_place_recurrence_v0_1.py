#!/usr/bin/env python3
"""Response-blind co-distribution census for the fresh 1,000-species TTF universe.

Reads only frozen species metadata and edge-midpoint geometry. It never reads
sequence lines, nucleotide identity, pairwise genetic distances or TTF outcomes.
"""
from __future__ import annotations
import argparse, csv, hashlib, json
from pathlib import Path
import numpy as np
from scipy.spatial import cKDTree

ARCHIVE_SHA="5a0fd9ac25893c749d14186fbcce4a46b99163c9d810b36e40eebce7bece61a5"
CSV_SHA="08dca43252ea93d43f3095f17b0267ec7ff049f0ec8b3e9c412dd784319ab575"
NPZ_SHA="d0776b4b5646e999f0c09ad937e4050fafe6fec018f17f0d37cd2c5b5f12a79b"
RADIUS=500.0
MIN_SYM=0.50
MIN_SOURCES=5

def sha256_path(path: Path) -> str:
    h=hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda:f.read(1<<20), b""):
            h.update(block)
    return h.hexdigest()

def hkey(tag: str, species: str) -> str:
    return hashlib.sha256(f"{tag}|{ARCHIVE_SHA}|{species}".encode()).hexdigest()

def digest_labels(labels) -> str:
    payload="".join(f"{x}\n" for x in sorted(map(str,labels)))
    return hashlib.sha256(payload.encode()).hexdigest()

def digest_pairs(pairs) -> str:
    payload="".join(f"{a}\t{b}\n" for a,b in sorted(pairs))
    return hashlib.sha256(payload.encode()).hexdigest()

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--candidate-csv", type=Path, required=True)
    ap.add_argument("--compact-npz", type=Path, required=True)
    ap.add_argument("--output", type=Path, required=True)
    args=ap.parse_args()
    if sha256_path(args.candidate_csv)!=CSV_SHA: raise SystemExit("candidate CSV hash mismatch")
    if sha256_path(args.compact_npz)!=NPZ_SHA: raise SystemExit("compact NPZ hash mismatch")

    with args.candidate_csv.open(newline="", encoding="utf-8") as f:
        rows=list(csv.DictReader(f))
    meta={r["species"]:r for r in rows}

    z=np.load(args.compact_npz, allow_pickle=False)
    names=np.asarray(z["species_order"]).astype(str)
    offsets=np.asarray(z["edge_offsets"], dtype=np.int64)
    points=np.asarray(z["edge_midpoints_ecef_km"], dtype=float)
    n_edges=np.asarray(z["n_edges"], dtype=np.int64)
    if len(names)!=1000 or offsets[-1]!=len(points): raise SystemExit("geometry shape drift")
    if set(names)!=set(meta): raise SystemExit("candidate/geometry species drift")

    centers=np.empty((len(names),3)); radii=np.empty(len(names))
    for i in range(len(names)):
        p=points[offsets[i]:offsets[i+1]]
        c=p.mean(axis=0); centers[i]=c
        radii[i]=np.sqrt(np.sum((p-c)**2,axis=1)).max()
    delta=centers[:,None,:]-centers[None,:,:]
    center_dist=np.sqrt(np.sum(delta*delta,axis=2))
    candidate=center_dist <= radii[:,None]+radii[None,:]+RADIUS
    np.fill_diagonal(candidate, False)

    coverage=np.zeros((len(names),len(names)), dtype=np.float32)  # target, source
    trees=[None]*len(names)
    for source in range(len(names)):
        targets=np.flatnonzero(candidate[:,source])
        if not len(targets): continue
        query=np.concatenate([points[offsets[t]:offsets[t+1]] for t in targets], axis=0)
        tree=cKDTree(points[offsets[source]:offsets[source+1]])
        dist,_=tree.query(query,k=1,distance_upper_bound=RADIUS,workers=1)
        lengths=n_edges[targets]
        starts=np.r_[0,np.cumsum(lengths)[:-1]]
        hits=np.add.reduceat((dist<=RADIUS).astype(np.int32), starts)
        coverage[targets,source]=hits/lengths

    symmetric=np.minimum(coverage,coverage.T)
    np.fill_diagonal(symmetric,0.0)
    codistributed=symmetric>=MIN_SYM

    order=sorted(range(len(names)), key=lambda i:(hkey("codistributed-place-v0.1-panel",names[i]),names[i]))
    dev=np.asarray(order[:500]); conf=np.asarray(order[500:])

    def split(panel, tag):
        x=sorted(panel.tolist(), key=lambda i:(hkey(tag,names[i]),names[i]))
        return np.asarray(x[:250]),np.asarray(x[250:])

    dtrain,deval=split(dev,"codistributed-place-v0.1-development-split")
    ctrain,ceval=split(conf,"codistributed-place-v0.1-confirmatory-split")

    def summarize(panel,train,evals):
        counts=codistributed[np.ix_(evals,train)].sum(axis=1).astype(int)
        supported=evals[counts>=MIN_SOURCES]
        pairs=[(names[s],names[t]) for t in supported for s in train if codistributed[t,s]]
        order_counts={}
        for t in supported:
            key=meta[names[t]]["order"]
            order_counts[key]=order_counts.get(key,0)+1
        q=np.quantile(counts,[0,.1,.25,.5,.75,.9,1])
        return {
          "panel_species":int(len(panel)),"train_species":int(len(train)),"eval_species":int(len(evals)),
          "panel_digest_sha256":digest_labels(names[panel]),
          "train_digest_sha256":digest_labels(names[train]),
          "eval_digest_sha256":digest_labels(names[evals]),
          "within_panel_codistributed_pairs":int(np.count_nonzero(np.triu(codistributed[np.ix_(panel,panel)],1))),
          "eval_targets_ge5_sources":int(np.count_nonzero(counts>=MIN_SOURCES)),
          "eval_targets_ge10_sources":int(np.count_nonzero(counts>=10)),
          "eval_source_count_quantiles":dict(zip(["min","q10","q25","median","q75","q90","max"],map(float,q))),
          "supported_eval_species_digest_sha256":digest_labels(names[supported]),
          "supported_dyads":len(pairs),
          "supported_dyad_digest_sha256":digest_pairs(pairs),
          "supported_eval_order_counts":dict(sorted(order_counts.items())),
        }

    all_pairs=[(names[i],names[j]) for i in range(len(names)) for j in range(i+1,len(names)) if codistributed[i,j]]
    out={
      "schema":"ttf_genetic_codistributed_place_response_blind_census_v0.1",
      "inputs":{"candidate_csv_sha256":CSV_SHA,"compact_geometry_npz_sha256":NPZ_SHA,"species":1000,"edge_midpoints":int(len(points))},
      "rule":{"radius_km":RADIUS,"symmetric_minimum_coverage":MIN_SYM,"minimum_sources":MIN_SOURCES},
      "all_1000":{"qualifying_undirected_pairs":len(all_pairs),"pair_graph_digest_sha256":digest_pairs(all_pairs)},
      "development":summarize(dev,dtrain,deval),
      "confirmatory":summarize(conf,ctrain,ceval),
      "outcome_firewall":{"nucleotide_identity_opened":False,"pairwise_genetic_distances_opened":False,"empirical_place_recurrence_statistic_opened":False},
    }
    args.output.write_text(json.dumps(out,indent=2,sort_keys=True)+"\n",encoding="utf-8")

if __name__=="__main__":
    main()
