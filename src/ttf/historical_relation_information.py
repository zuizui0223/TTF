from __future__ import annotations

from dataclasses import asdict, dataclass

import numpy as np

from .relational_qualification import (
    dyadic_signal_support,
    residualize_primary_relation,
)


@dataclass(frozen=True)
class HistoricalRelationEcology:
    dyads: int
    spearman_hist_vs_current: float
    current_q75: float
    historical_q25: float
    historical_q75: float
    high_current_dyads: int
    high_current_low_history_dyads: int
    high_current_high_history_dyads: int
    high_current_low_history_fraction: float
    high_current_high_history_fraction: float


def _zscore(values: np.ndarray) -> np.ndarray:
    x=np.asarray(values,float)
    if x.ndim!=1 or len(x)==0 or not np.isfinite(x).all():
        raise ValueError("expected finite non-empty vector")
    sd=float(x.std())
    if sd<=np.finfo(float).eps:
        raise ValueError("constant vector")
    return (x-float(x.mean()))/sd


def _rank(values: np.ndarray) -> np.ndarray:
    x=np.asarray(values,float)
    order=np.argsort(x,kind="stable")
    out=np.empty(len(x),float)
    sx=x[order]
    starts=np.r_[0,1+np.flatnonzero(sx[1:]!=sx[:-1])]
    stops=np.r_[starts[1:],len(x)]
    for start,stop in zip(starts,stops):
        out[order[start:stop]]=0.5*((start+1)+stop)
    return out


def _spearman(a: np.ndarray,b: np.ndarray) -> float:
    x=_rank(a); y=_rank(b)
    x-=x.mean(); y-=y.mean()
    den=float(np.sqrt(np.dot(x,x)*np.dot(y,y)))
    return 0.0 if den<=np.finfo(float).eps else float(np.dot(x,y)/den)


def historical_relation_information(
    *,
    species_order: np.ndarray,
    class_name: np.ndarray,
    order_name: np.ndarray,
    family_name: np.ndarray,
    occurrence_n: np.ndarray,
    source_index: np.ndarray,
    target_index: np.ndarray,
    r_hist: np.ndarray,
    r_current: np.ndarray,
) -> dict[str,object]:
    """Characterize historical-climate relation information without a response.

    This intentionally stops before geographic-opportunity and genetic-response
    stages. It asks how much historical similarity remains after present climate,
    lineage, locality-count imbalance, and exact source/target identity.
    """
    species=np.asarray(species_order).astype(str)
    cls=np.asarray(class_name).astype(str)
    order=np.asarray(order_name).astype(str)
    family=np.asarray(family_name).astype(str)
    occ=np.asarray(occurrence_n,float)
    s=np.asarray(source_index,np.int64)
    t=np.asarray(target_index,np.int64)
    rh=np.asarray(r_hist,float)
    rc=np.asarray(r_current,float)

    n=len(species)
    if not (
        cls.shape==(n,) and order.shape==(n,) and family.shape==(n,)
        and occ.shape==(n,)
    ):
        raise ValueError("species metadata arrays must align")
    if not (
        s.ndim==t.ndim==rh.ndim==rc.ndim==1
        and len(s)==len(t)==len(rh)==len(rc)
        and len(s)>0
    ):
        raise ValueError("dyad arrays must align and be non-empty")
    if np.any(s<0) or np.any(t<0) or np.any(s>=n) or np.any(t>=n):
        raise ValueError("dyad index outside species_order")
    if not (np.isfinite(occ).all() and np.isfinite(rh).all() and np.isfinite(rc).all()):
        raise ValueError("all numeric inputs must be finite")
    if np.any(occ<=0):
        raise ValueError("occurrence counts must be positive")

    same_class=(cls[s]==cls[t]).astype(float)
    same_order=(order[s]==order[t]).astype(float)
    same_family=(family[s]==family[t]).astype(float)
    locality_ratio=occ[s]/occ[t]

    primary=_zscore(rh)
    controls=np.column_stack((
        _zscore(rc),
        same_class-same_class.mean(),
        same_order-same_order.mean(),
        same_family-same_family.mean(),
        _zscore(np.abs(np.log(locality_ratio))),
    ))
    residual,nonredundancy=residualize_primary_relation(
        s,t,primary,controls
    )
    support=dyadic_signal_support(s,t,residual)

    current_q75=float(np.quantile(rc,0.75))
    hist_q25=float(np.quantile(rh,0.25))
    hist_q75=float(np.quantile(rh,0.75))
    high_current=rc>=current_q75
    high_current_n=int(np.count_nonzero(high_current))
    low_history=high_current & (rh<=hist_q25)
    high_history=high_current & (rh>=hist_q75)
    denom=max(high_current_n,1)

    ecology=HistoricalRelationEcology(
        dyads=int(len(rh)),
        spearman_hist_vs_current=_spearman(rh,rc),
        current_q75=current_q75,
        historical_q25=hist_q25,
        historical_q75=hist_q75,
        high_current_dyads=high_current_n,
        high_current_low_history_dyads=int(np.count_nonzero(low_history)),
        high_current_high_history_dyads=int(np.count_nonzero(high_history)),
        high_current_low_history_fraction=float(np.count_nonzero(low_history)/denom),
        high_current_high_history_fraction=float(np.count_nonzero(high_history)/denom),
    )
    return {
        "ecological_geometry":asdict(ecology),
        "nonredundancy":asdict(nonredundancy),
        "dyadic_signal_support":asdict(support),
        "conditioning_set":[
            "current_climate_similarity",
            "same_class",
            "same_order",
            "same_family",
            "absolute_log_occurrence_count_ratio",
            "source_fixed_effect",
            "target_fixed_effect",
        ],
        "excluded_at_this_stage":[
            "geographic_opportunity",
            "genetic_response",
        ],
    }
