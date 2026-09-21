from __future__ import annotations

from dataclasses import dataclass
import hashlib
from typing import Mapping, Sequence

import numpy as np


@dataclass(frozen=True)
class SpeciesTransferabilityDecomposition:
    source_names: tuple[str, ...]
    target_names: tuple[str, ...]
    source_effect: np.ndarray
    target_effect: np.ndarray
    source_se: np.ndarray
    target_se: np.ndarray
    source_variance: float
    target_variance: float
    residual_variance: float
    source_fraction: float
    target_fraction: float
    residual_fraction: float


@dataclass(frozen=True)
class SpeciesTransferabilityCV:
    folds: int
    geometry_only_rmse: float
    species_property_rmse: float
    rmse_improvement: float
    geometry_only_correlation: float
    species_property_correlation: float


def _standardize_columns(x: np.ndarray) -> np.ndarray:
    x=np.asarray(x,float)
    if x.ndim!=2:
        raise ValueError("x must be a matrix")
    center=np.mean(x,axis=0)
    scale=np.std(x,axis=0,ddof=0)
    scale[scale<=np.sqrt(np.finfo(float).eps)]=1.0
    return (x-center)/scale


def geometry_residual(
    pair_transfer: Sequence[float] | np.ndarray,
    geometry_covariates: np.ndarray,
) -> tuple[np.ndarray,np.ndarray]:
    y=np.asarray(pair_transfer,float)
    x=_standardize_columns(np.asarray(geometry_covariates,float))
    if y.ndim!=1 or len(y)!=len(x) or len(y)<8:
        raise ValueError("pair response/design drift")
    if not np.isfinite(y).all() or not np.isfinite(x).all():
        raise ValueError("non-finite pair data")
    X=np.column_stack([np.ones(len(y)),x])
    beta=np.linalg.lstsq(X,y,rcond=None)[0]
    residual=y-X@beta
    return residual,beta


def _same_label_covariance(values: np.ndarray, labels: np.ndarray) -> float:
    total=0.0
    pairs=0
    for label in np.unique(labels):
        v=np.asarray(values[labels==label],float)
        if len(v)<2:
            continue
        # sum over unordered products without O(n^2) materialization
        total += 0.5*(float(np.sum(v))**2-float(np.dot(v,v)))
        pairs += len(v)*(len(v)-1)//2
    return 0.0 if pairs==0 else float(total/pairs)


def moment_variance_components(
    residual: Sequence[float] | np.ndarray,
    source_index: Sequence[int] | np.ndarray,
    target_index: Sequence[int] | np.ndarray,
) -> tuple[float,float,float]:
    r=np.asarray(residual,float)
    s=np.asarray(source_index)
    t=np.asarray(target_index)
    if r.ndim!=1 or len(r)!=len(s) or len(r)!=len(t):
        raise ValueError("pair arrays drift")
    centered=r-float(np.mean(r))
    source=max(0.0,_same_label_covariance(centered,s))
    target=max(0.0,_same_label_covariance(centered,t))
    total=float(np.mean(centered*centered))
    error=max(np.finfo(float).eps,total-source-target)
    return source,target,error


def _blup(
    residual: np.ndarray,
    source_index: np.ndarray,
    target_index: np.ndarray,
    n_source: int,
    n_target: int,
    source_variance: float,
    target_variance: float,
    residual_variance: float,
    *,
    compute_se: bool = True,
) -> tuple[np.ndarray,np.ndarray,np.ndarray,np.ndarray]:
    n=len(residual)
    q=n_source+n_target
    Z=np.zeros((n,q),float)
    Z[np.arange(n),source_index]=1.0
    Z[np.arange(n),n_source+target_index]=1.0
    precision=(Z.T@Z)/residual_variance
    # A tiny floor represents an effectively zero random-effect variance.
    sv=max(source_variance,1e-12)
    tv=max(target_variance,1e-12)
    precision[np.arange(n_source),np.arange(n_source)] += 1.0/sv
    j=np.arange(n_target)+n_source
    precision[j,j] += 1.0/tv
    rhs=(Z.T@residual)/residual_variance
    effects=np.linalg.solve(precision,rhs)
    if compute_se:
        covariance=np.linalg.inv(precision)
        se=np.sqrt(np.maximum(np.diag(covariance),0.0))
    else:
        se=np.full(q,np.nan,dtype=float)
    return effects[:n_source],effects[n_source:],se[:n_source],se[n_source:]


def decompose_species_transferability(
    pair_transfer: Sequence[float] | np.ndarray,
    geometry_covariates: np.ndarray,
    source_species: Sequence[str],
    target_species: Sequence[str],
) -> SpeciesTransferabilityDecomposition:
    y=np.asarray(pair_transfer,float)
    source=np.asarray(list(map(str,source_species)),dtype=object)
    target=np.asarray(list(map(str,target_species)),dtype=object)
    if len(y)!=len(source) or len(y)!=len(target):
        raise ValueError("pair labels drift")
    source_names=tuple(sorted(set(source)))
    target_names=tuple(sorted(set(target)))
    smap={name:i for i,name in enumerate(source_names)}
    tmap={name:i for i,name in enumerate(target_names)}
    si=np.asarray([smap[x] for x in source],np.int64)
    ti=np.asarray([tmap[x] for x in target],np.int64)
    residual,_=geometry_residual(y,geometry_covariates)
    sv,tv,ev=moment_variance_components(residual,si,ti)
    us,vt,use,tse=_blup(residual,si,ti,len(source_names),len(target_names),sv,tv,ev)
    total=sv+tv+ev
    return SpeciesTransferabilityDecomposition(
        source_names=source_names,target_names=target_names,
        source_effect=us,target_effect=vt,source_se=use,target_se=tse,
        source_variance=sv,target_variance=tv,residual_variance=ev,
        source_fraction=float(sv/total),target_fraction=float(tv/total),
        residual_fraction=float(ev/total),
    )


def deterministic_pair_folds(
    source_species: Sequence[str],
    target_species: Sequence[str],
    *,
    folds: int=10,
    tag: str="lepidoptera-species-transferability-v0.1",
) -> np.ndarray:
    if folds<2:
        raise ValueError("folds must be >=2")
    out=[]
    for source,target in zip(source_species,target_species):
        digest=hashlib.sha256(f"{tag}|{source}|{target}".encode()).digest()
        out.append(int.from_bytes(digest[:8],"big")%int(folds))
    return np.asarray(out,np.int64)


def _corr(a: np.ndarray,b: np.ndarray) -> float:
    aa=a-float(np.mean(a)); bb=b-float(np.mean(b))
    denom=float(np.sqrt(np.dot(aa,aa)*np.dot(bb,bb)))
    return 0.0 if denom<=np.finfo(float).tiny else float(np.dot(aa,bb)/denom)


def crossvalidated_species_property_gain(
    pair_transfer: Sequence[float] | np.ndarray,
    geometry_covariates: np.ndarray,
    source_species: Sequence[str],
    target_species: Sequence[str],
    *,
    folds: int=10,
) -> SpeciesTransferabilityCV:
    y=np.asarray(pair_transfer,float)
    g=np.asarray(geometry_covariates,float)
    source=np.asarray(list(map(str,source_species)),dtype=object)
    target=np.asarray(list(map(str,target_species)),dtype=object)
    fold=deterministic_pair_folds(source,target,folds=folds)
    pred_geo=np.empty(len(y),float)
    pred_species=np.empty(len(y),float)
    for k in range(int(folds)):
        train=fold!=k; test=~train
        xtrain=_standardize_columns(g[train])
        center=np.mean(g[train],axis=0)
        scale=np.std(g[train],axis=0,ddof=0)
        scale[scale<=np.sqrt(np.finfo(float).eps)]=1.0
        xtest=(g[test]-center)/scale
        Xtrain=np.column_stack([np.ones(np.count_nonzero(train)),xtrain])
        Xtest=np.column_stack([np.ones(np.count_nonzero(test)),xtest])
        beta=np.linalg.lstsq(Xtrain,y[train],rcond=None)[0]
        pred_geo[test]=Xtest@beta
        residual_train=y[train]-Xtrain@beta

        source_names=tuple(sorted(set(source[train])))
        target_names=tuple(sorted(set(target[train])))
        smap={name:i for i,name in enumerate(source_names)}
        tmap={name:i for i,name in enumerate(target_names)}
        si=np.asarray([smap[x] for x in source[train]],np.int64)
        ti=np.asarray([tmap[x] for x in target[train]],np.int64)
        sv,tv,ev=moment_variance_components(residual_train,si,ti)
        us,vt,_,_=_blup(
            residual_train,si,ti,len(source_names),len(target_names),sv,tv,ev,
            compute_se=False,
        )
        extra=np.zeros(np.count_nonzero(test),float)
        test_sources=source[test]; test_targets=target[test]
        for i,(sname,tname) in enumerate(zip(test_sources,test_targets)):
            if sname in smap:
                extra[i]+=us[smap[sname]]
            if tname in tmap:
                extra[i]+=vt[tmap[tname]]
        pred_species[test]=pred_geo[test]+extra
    rmse_geo=float(np.sqrt(np.mean((y-pred_geo)**2)))
    rmse_species=float(np.sqrt(np.mean((y-pred_species)**2)))
    return SpeciesTransferabilityCV(
        folds=int(folds),
        geometry_only_rmse=rmse_geo,
        species_property_rmse=rmse_species,
        rmse_improvement=float(rmse_geo-rmse_species),
        geometry_only_correlation=_corr(y,pred_geo),
        species_property_correlation=_corr(y,pred_species),
    )


__all__=[
    "SpeciesTransferabilityCV",
    "SpeciesTransferabilityDecomposition",
    "crossvalidated_species_property_gain",
    "decompose_species_transferability",
    "deterministic_pair_folds",
    "geometry_residual",
    "moment_variance_components",
]
