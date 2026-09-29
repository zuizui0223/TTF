#!/usr/bin/env python3
from __future__ import annotations

import argparse
import hashlib
from pathlib import Path

import numpy as np

from ttf.palearctic_lgm import (
    fit_weighted_ridge_logistic,
    predict_relative_suitability,
    schoener_d,
)


def seed_for(replicate: int, species: str) -> int:
    token=f"palearctic-insect-lgm-bootstrap-v0.1|{replicate}|{species}"
    return int(hashlib.sha256(token.encode()).hexdigest()[:16],16) % (2**63-1)


def main() -> int:
    ap=argparse.ArgumentParser()
    ap.add_argument("--design",type=Path,required=True)
    ap.add_argument("--replicate",type=int,required=True)
    ap.add_argument("--output",type=Path,required=True)
    args=ap.parse_args()
    if not (0 <= args.replicate < 100):
        raise RuntimeError("bootstrap replicate must be in [0,99]")

    data=np.load(args.design,allow_pickle=False)
    species=list(map(str,data["species_order"]))
    offsets=np.asarray(data["presence_offsets"],dtype=np.int64)
    pres=np.asarray(data["presence_features"],dtype=float)
    bg=np.asarray(data["background_features"],dtype=float)
    lgm_grid=np.asarray(data["lgm_grid_features"],dtype=float)
    source=np.asarray(data["source_index"],dtype=np.int64)
    target=np.asarray(data["target_index"],dtype=np.int64)

    surfaces=[]
    for i,name in enumerate(species):
        x=pres[offsets[i]:offsets[i+1]]
        rng=np.random.default_rng(seed_for(args.replicate,name))
        resampled=x[rng.integers(0,len(x),size=len(x))]
        fit=fit_weighted_ridge_logistic(
            resampled,bg,
            ridge_lambda=1.0,
            presence_total_weight=.5,
            background_total_weight=.5,
            max_iter=2000,
        )
        if not fit.converged:
            raise RuntimeError(
                f"bootstrap SDM failed to converge: replicate={args.replicate}, species={name}"
            )
        surfaces.append(predict_relative_suitability(lgm_grid,fit))
    s=np.vstack(surfaces)
    relation=np.asarray([
        schoener_d(s[int(i)],s[int(j)])
        for i,j in zip(source,target)
    ],dtype=float)
    args.output.parent.mkdir(parents=True,exist_ok=True)
    np.save(args.output,relation,allow_pickle=False)
    print({
        "replicate":args.replicate,
        "species":len(species),
        "dyads":len(relation),
        "relation_mean":float(relation.mean()),
    })
    return 0


if __name__=="__main__":
    raise SystemExit(main())
