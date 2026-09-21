from __future__ import annotations

from dataclasses import dataclass
from typing import Mapping, Sequence
import numpy as np


GEOMETRY_COVARIATES = (
    "coverage",
    "centroid_distance",
    "edge_count_ratio",
    "locality_count_ratio",
)


@dataclass(frozen=True)
class TraitGradientPairDesign:
    target: np.ndarray
    source: np.ndarray
    similarity: np.ndarray
    residual_similarity: np.ndarray
    target_names: tuple[str, ...]
    source_names: tuple[str, ...]


def _standardize(x: np.ndarray) -> np.ndarray:
    x=np.asarray(x,float)
    sd=float(np.std(x))
    return np.zeros_like(x) if sd == 0 else (x-float(np.mean(x)))/sd


def continuous_similarity(values: Mapping[str,float], a: str, b: str) -> float:
    names=tuple(sorted(values))
    z=dict(zip(names,_standardize(np.asarray([np.log1p(max(0.0,float(values[n]))) for n in names]))))
    return float(np.exp(-abs(z[a]-z[b])))


def categorical_similarity(values: Mapping[str,str], a: str, b: str) -> float:
    return float(str(values[a]).strip() == str(values[b]).strip())


def habitat_similarity(values: Mapping[str,Sequence[str]], a: str, b: str) -> float:
    xa=tuple(map(str,values[a])); xb=tuple(map(str,values[b]))
    if len(xa)!=4 or len(xb)!=4: raise ValueError("habitat requires four frozen axes")
    return float(np.mean([u.strip()==v.strip() for u,v in zip(xa,xb)]))


def composite_trait_similarity(
    target: str,
    source: str,
    *,
    wing_size: Mapping[str,float],
    voltinism: Mapping[str,str],
    host_breadth: Mapping[str,float],
    habitat: Mapping[str,Sequence[str]],
) -> float:
    return float(np.mean([
        continuous_similarity(wing_size,target,source),
        categorical_similarity(voltinism,target,source),
        continuous_similarity(host_breadth,target,source),
        habitat_similarity(habitat,target,source),
    ]))


def residualize_within_target(similarity: np.ndarray, geometry: np.ndarray, target: np.ndarray) -> np.ndarray:
    y=np.asarray(similarity,float); x=np.asarray(geometry,float); t=np.asarray(target)
    if x.ndim!=2 or x.shape[1]!=4 or len(y)!=len(x) or len(t)!=len(y): raise ValueError("pair arrays drift")
    out=np.empty_like(y)
    for label in np.unique(t):
        idx=np.flatnonzero(t==label)
        X=np.column_stack([np.ones(len(idx)),x[idx]])
        beta=np.linalg.lstsq(X,y[idx],rcond=None)[0]
        out[idx]=y[idx]-X@beta
    return out


def equal_target_gradient(
    transfer_contribution: np.ndarray,
    residual_similarity: np.ndarray,
    target: np.ndarray,
) -> tuple[float,dict[str,float]]:
    y=np.asarray(transfer_contribution,float); x=np.asarray(residual_similarity,float); t=np.asarray(target)
    slopes={}
    for label in np.unique(t):
        idx=np.flatnonzero(t==label); xx=x[idx]; yy=y[idx]
        denom=float(np.dot(xx,xx))
        if len(idx)<3 or denom<=0: continue
        slopes[str(label)]=float(np.dot(xx,yy)/denom)
    if not slopes: raise ValueError("no target-specific gradient is estimable")
    return float(np.mean(list(slopes.values()))),slopes
