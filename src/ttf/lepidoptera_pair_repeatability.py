from __future__ import annotations

from dataclasses import dataclass
import hashlib
from typing import Mapping

import numpy as np


SPLIT_TAG="pair-repeatability-v01"


@dataclass(frozen=True)
class PreparedPairRepeatability:
    pair_index: np.ndarray
    target_index: np.ndarray
    source_index: np.ndarray
    target_names: tuple[str,...]
    source_names: tuple[str,...]
    pair_covariates: np.ndarray
    residualized_covariates: np.ndarray
    source_group: np.ndarray
    target_group: np.ndarray
    source_group_count: int
    target_group_count: int
    half_a_pair_id: np.ndarray
    half_b_pair_id: np.ndarray
    half_a_global_target: np.ndarray
    half_a_global_source: np.ndarray
    half_b_global_target: np.ndarray
    half_b_global_source: np.ndarray
    half_a_rows: np.ndarray
    half_b_rows: np.ndarray


@dataclass(frozen=True)
class PairRepeatabilityScore:
    statistic: float
    covariance: float
    half_a_variance: float
    half_b_variance: float
    pair_variance_fraction: float
    retained_pairs: int


def _row_hash(target_name:str,source_name:str,target_edge:int,source_edge:int)->bytes:
    payload=(
        f"{SPLIT_TAG}|{target_name}|{source_name}|"
        f"{int(target_edge)}|{int(source_edge)}"
    ).encode("utf-8")
    return hashlib.sha256(payload).digest()


def _two_way_demean(
    values: np.ndarray,
    source_group: np.ndarray,
    target_group: np.ndarray,
    n_source: int,
    n_target: int,
    *,
    max_iter: int=200,
    tolerance: float=1e-12,
) -> np.ndarray:
    x=np.asarray(values,float).copy()
    vector=x.ndim==1
    if vector:
        x=x[:,None]
    if x.ndim!=2 or len(x)!=len(source_group) or len(x)!=len(target_group):
        raise ValueError("two-way demeaning arrays drift")
    sc=np.bincount(source_group,minlength=n_source).astype(float)
    tc=np.bincount(target_group,minlength=n_target).astype(float)
    if np.any(sc<=0) or np.any(tc<=0):
        raise ValueError("empty source/target group")
    for _ in range(int(max_iter)):
        before=x.copy()
        for j in range(x.shape[1]):
            sm=np.bincount(source_group,weights=x[:,j],minlength=n_source)/sc
            x[:,j]-=sm[source_group]
            tm=np.bincount(target_group,weights=x[:,j],minlength=n_target)/tc
            x[:,j]-=tm[target_group]
        if float(np.max(np.abs(x-before)))<=float(tolerance):
            break
    return x[:,0] if vector else x


def _residualize_half(
    y: np.ndarray,
    prepared: PreparedPairRepeatability,
) -> np.ndarray:
    yd=_two_way_demean(
        np.asarray(y,float),
        prepared.source_group,
        prepared.target_group,
        prepared.source_group_count,
        prepared.target_group_count,
    )
    x=prepared.residualized_covariates
    if x.shape[1]:
        beta=np.linalg.lstsq(x,yd,rcond=None)[0]
        yd=yd-x@beta
    yd=_two_way_demean(
        yd,
        prepared.source_group,
        prepared.target_group,
        prepared.source_group_count,
        prepared.target_group_count,
    )
    return yd


def _pair_correlation(
    flat_response: np.ndarray,
    pair_id: np.ndarray,
    global_target: np.ndarray,
    global_source: np.ndarray,
    n_pairs: int,
) -> np.ndarray:
    x=flat_response[global_target]
    y=flat_response[global_source]
    n=np.bincount(pair_id,minlength=n_pairs).astype(float)
    sx=np.bincount(pair_id,weights=x,minlength=n_pairs)
    sy=np.bincount(pair_id,weights=y,minlength=n_pairs)
    sxx=np.bincount(pair_id,weights=x*x,minlength=n_pairs)
    syy=np.bincount(pair_id,weights=y*y,minlength=n_pairs)
    sxy=np.bincount(pair_id,weights=x*y,minlength=n_pairs)
    cov=sxy-sx*sy/n
    den=np.sqrt(np.maximum((sxx-sx*sx/n)*(syy-sy*sy/n),0.0))
    out=np.zeros(n_pairs,float)
    ok=den>np.sqrt(np.finfo(float).eps)
    out[ok]=cov[ok]/den[ok]
    return np.clip(out,-0.999999,0.999999)


def prepare_pair_repeatability(
    design_npz,
    edge_counts: np.ndarray,
    *,
    minimum_rows_per_half: int=8,
) -> PreparedPairRepeatability:
    names=tuple(map(str,design_npz["species_order"]))
    target_all=np.asarray(design_npz["target_index"],np.int64)
    source_all=np.asarray(design_npz["source_index"],np.int64)
    lengths=np.asarray(design_npz["alignment_length"],np.int64)
    offsets=np.asarray(design_npz["alignment_offset"],np.int64)
    at=np.asarray(design_npz["alignment_target"],np.int64)
    ass=np.asarray(design_npz["alignment_source"],np.int64)
    expected=np.r_[0,np.cumsum(lengths)[:-1]]
    if not np.array_equal(offsets,expected):
        raise RuntimeError("alignment offset drift")
    keep=np.flatnonzero(lengths>=2*int(minimum_rows_per_half))
    if len(keep)<100:
        raise RuntimeError("insufficient repeatability pairs")

    source_values=tuple(sorted(set(map(int,source_all[keep]))))
    target_values=tuple(sorted(set(map(int,target_all[keep]))))
    smap={v:i for i,v in enumerate(source_values)}
    tmap={v:i for i,v in enumerate(target_values)}
    source_group=np.asarray([smap[int(v)] for v in source_all[keep]],np.int64)
    target_group=np.asarray([tmap[int(v)] for v in target_all[keep]],np.int64)

    nuisance=np.asarray(design_npz["nuisance_covariates"],float)[keep]
    raw=np.asarray(design_npz["raw_host_resource_jaccard"],float)[keep,None]
    log_rows=np.log1p(lengths[keep].astype(float))[:,None]
    covariates=np.column_stack([nuisance,raw,log_rows])
    covariates=(covariates-np.mean(covariates,axis=0))/np.where(
        np.std(covariates,axis=0)>np.sqrt(np.finfo(float).eps),
        np.std(covariates,axis=0),
        1.0,
    )
    covariates=_two_way_demean(
        covariates,source_group,target_group,len(source_values),len(target_values)
    )

    edge_counts=np.asarray(edge_counts,np.int64)
    if len(edge_counts)!=len(names):
        raise ValueError("edge_counts drift")
    edge_offsets=np.r_[0,np.cumsum(edge_counts)]

    a_pid=[]; b_pid=[]; a_gt=[]; a_gs=[]; b_gt=[]; b_gs=[]
    a_rows=[]; b_rows=[]
    for new_id,old_id in enumerate(keep):
        start=int(offsets[old_id]); n=int(lengths[old_id])
        te=at[start:start+n]; se=ass[start:start+n]
        tname=names[int(target_all[old_id])]
        sname=names[int(source_all[old_id])]
        order=sorted(
            range(n),
            key=lambda j:_row_hash(tname,sname,int(te[j]),int(se[j])),
        )
        cut=n//2
        aa=np.asarray(order[:cut],np.int64)
        bb=np.asarray(order[cut:],np.int64)
        if len(aa)<minimum_rows_per_half or len(bb)<minimum_rows_per_half:
            raise RuntimeError("split-half support drift")
        a_pid.extend([new_id]*len(aa)); b_pid.extend([new_id]*len(bb))
        a_gt.extend((edge_offsets[int(target_all[old_id])]+te[aa]).tolist())
        a_gs.extend((edge_offsets[int(source_all[old_id])]+se[aa]).tolist())
        b_gt.extend((edge_offsets[int(target_all[old_id])]+te[bb]).tolist())
        b_gs.extend((edge_offsets[int(source_all[old_id])]+se[bb]).tolist())
        a_rows.append(len(aa)); b_rows.append(len(bb))

    return PreparedPairRepeatability(
        pair_index=keep,
        target_index=target_all[keep],
        source_index=source_all[keep],
        target_names=tuple(names[int(i)] for i in target_all[keep]),
        source_names=tuple(names[int(i)] for i in source_all[keep]),
        pair_covariates=np.column_stack([
            np.asarray(design_npz["nuisance_covariates"],float)[keep],
            np.asarray(design_npz["raw_host_resource_jaccard"],float)[keep,None],
            log_rows,
        ]),
        residualized_covariates=covariates,
        source_group=source_group,
        target_group=target_group,
        source_group_count=len(source_values),
        target_group_count=len(target_values),
        half_a_pair_id=np.asarray(a_pid,np.int64),
        half_b_pair_id=np.asarray(b_pid,np.int64),
        half_a_global_target=np.asarray(a_gt,np.int64),
        half_a_global_source=np.asarray(a_gs,np.int64),
        half_b_global_target=np.asarray(b_gt,np.int64),
        half_b_global_source=np.asarray(b_gs,np.int64),
        half_a_rows=np.asarray(a_rows,np.int64),
        half_b_rows=np.asarray(b_rows,np.int64),
    )


def score_pair_repeatability(
    response_by_species: Mapping[str,np.ndarray],
    species_order: tuple[str,...],
    prepared: PreparedPairRepeatability,
) -> PairRepeatabilityScore:
    if set(response_by_species)!=set(species_order):
        raise ValueError("response species differ from repeatability design")
    flat=np.concatenate([
        np.asarray(response_by_species[name],float) for name in species_order
    ])
    n=len(prepared.pair_index)
    ca=_pair_correlation(
        flat,prepared.half_a_pair_id,
        prepared.half_a_global_target,prepared.half_a_global_source,n
    )
    cb=_pair_correlation(
        flat,prepared.half_b_pair_id,
        prepared.half_b_global_target,prepared.half_b_global_source,n
    )
    za=np.arctanh(ca); zb=np.arctanh(cb)
    ra=_residualize_half(za,prepared)
    rb=_residualize_half(zb,prepared)
    aa=ra-float(np.mean(ra)); bb=rb-float(np.mean(rb))
    va=float(np.mean(aa*aa)); vb=float(np.mean(bb*bb))
    covariance=float(np.mean(aa*bb))
    denom=float(np.sqrt(va*vb))
    statistic=0.0 if denom<=np.finfo(float).tiny else float(covariance/denom)
    mean_var=.5*(va+vb)
    fraction=0.0 if mean_var<=np.finfo(float).tiny else float(covariance/mean_var)
    return PairRepeatabilityScore(
        statistic=statistic,
        covariance=covariance,
        half_a_variance=va,
        half_b_variance=vb,
        pair_variance_fraction=fraction,
        retained_pairs=n,
    )


__all__=[
    "PairRepeatabilityScore",
    "PreparedPairRepeatability",
    "prepare_pair_repeatability",
    "score_pair_repeatability",
]
