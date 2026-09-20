from __future__ import annotations

from dataclasses import dataclass
from typing import Mapping
import numpy as np

from .genetic_geometry import GeneticSamplingGeometry
from .lepidoptera_trait_gradient_v02 import (
    envelope_pvalues,
    target_equal_mean_correlation,
    target_fe_geometry_residual,
)
from .phylogatr_compact_ibd import (
    crossfit_ibd_residuals_compact,
    prepare_compact_crossfit_ibd_design,
)


@dataclass(frozen=True)
class EmpiricalTraitGradientScore:
    statistic: float
    p_value: float
    positive: bool
    per_target: Mapping[int, float]


def post_ibd_edge_response(
    genetic_distance: np.ndarray,
    geometry: GeneticSamplingGeometry,
    *,
    min_training_edges: int = 5,
) -> np.ndarray:
    nodes=np.asarray(geometry.edge_nodes,dtype=np.int64)
    coords=np.asarray(geometry.coordinates,dtype=float)
    geographic=np.linalg.norm(
        coords[nodes[:,0]]-coords[nodes[:,1]],axis=1
    )
    design=prepare_compact_crossfit_ibd_design(
        geographic,nodes,min_training_edges=int(min_training_edges)
    )
    result=crossfit_ibd_residuals_compact(
        np.asarray(genetic_distance,dtype=float),design
    )
    response=np.asarray(result.residual_turnover,dtype=float)
    if response.shape!=(geometry.n_edges,) or not np.isfinite(response).all():
        raise RuntimeError("post-IBD response drift")
    return response


def frozen_pair_surface(design_npz) -> tuple[np.ndarray,np.ndarray,np.ndarray,np.ndarray,np.ndarray]:
    pairs=np.asarray(design_npz["pairs"],dtype=float)
    target=pairs[:,0].astype(np.int64)
    source=pairs[:,1].astype(np.int64)
    geometry=np.asarray(pairs[:,3:7],dtype=float).copy()
    geometry[:,1]=np.log1p(geometry[:,1]/500.0)
    geometry[:,2]=np.log(geometry[:,2])
    geometry[:,3]=np.log(geometry[:,3])
    kernel=np.asarray(design_npz["geometry_kernel"],dtype=float)
    geometry=np.column_stack([geometry,kernel[target,source]])
    residual=target_fe_geometry_residual(
        np.asarray(pairs[:,2],dtype=float),geometry,target
    )
    lengths=pairs[:,8].astype(np.int64)
    offsets=pairs[:,7].astype(np.int64)
    expected=np.r_[0,np.cumsum(lengths)[:-1]]
    if not np.array_equal(offsets,expected):
        raise RuntimeError("alignment offset drift")
    pair_id=np.repeat(np.arange(len(pairs),dtype=np.int64),lengths)
    return pairs,target,source,residual,pair_id


def pair_transfer_congruence(
    response_by_species: Mapping[str,np.ndarray],
    species_order: tuple[str,...],
    pairs: np.ndarray,
    pair_id: np.ndarray,
    alignment_target: np.ndarray,
    alignment_source: np.ndarray,
) -> np.ndarray:
    edge_counts=np.asarray(
        [len(np.asarray(response_by_species[name])) for name in species_order],
        dtype=np.int64,
    )
    edge_offsets=np.r_[0,np.cumsum(edge_counts)]
    target=pairs[:,0].astype(np.int64)
    source=pairs[:,1].astype(np.int64)
    flat=np.concatenate([
        np.asarray(response_by_species[name],dtype=float) for name in species_order
    ])
    global_target=edge_offsets[target][pair_id]+np.asarray(alignment_target,dtype=np.int64)
    global_source=edge_offsets[source][pair_id]+np.asarray(alignment_source,dtype=np.int64)
    x=flat[global_target]
    y=flat[global_source]
    n=np.bincount(pair_id,minlength=len(pairs)).astype(float)
    sx=np.bincount(pair_id,weights=x,minlength=len(pairs))
    sy=np.bincount(pair_id,weights=y,minlength=len(pairs))
    sxx=np.bincount(pair_id,weights=x*x,minlength=len(pairs))
    syy=np.bincount(pair_id,weights=y*y,minlength=len(pairs))
    sxy=np.bincount(pair_id,weights=x*y,minlength=len(pairs))
    covariance=sxy-sx*sy/n
    denominator=np.sqrt(np.maximum(
        (sxx-sx*sx/n)*(syy-sy*sy/n),0.0
    ))
    correlation=np.zeros(len(pairs),dtype=float)
    valid=denominator>np.sqrt(np.finfo(float).eps)
    correlation[valid]=covariance[valid]/denominator[valid]
    return correlation


def score_empirical_trait_gradient(
    response_by_species: Mapping[str,np.ndarray],
    design_npz,
    *,
    private_reference: np.ndarray,
    geometry_reference: np.ndarray,
    alpha: float = 0.05,
) -> EmpiricalTraitGradientScore:
    names=tuple(map(str,design_npz["species_order"]))
    if set(response_by_species)!=set(names):
        raise ValueError("response species differ from survivor design")
    pairs,target,source,residual,pair_id=frozen_pair_surface(design_npz)
    pair_score=pair_transfer_congruence(
        response_by_species,names,pairs,pair_id,
        np.asarray(design_npz["alignment_target"],dtype=np.int64),
        np.asarray(design_npz["alignment_source"],dtype=np.int64),
    )
    statistic,per_target=target_equal_mean_correlation(
        pair_score,residual,target
    )
    p=float(envelope_pvalues(
        np.asarray([statistic],dtype=float),
        np.asarray(private_reference,dtype=float),
        np.asarray(geometry_reference,dtype=float),
    )[0])
    return EmpiricalTraitGradientScore(
        statistic=float(statistic),
        p_value=p,
        positive=bool(float(statistic)>0.0 and p<=float(alpha)),
        per_target=per_target,
    )


__all__=[
    "EmpiricalTraitGradientScore",
    "frozen_pair_surface",
    "pair_transfer_congruence",
    "post_ibd_edge_response",
    "score_empirical_trait_gradient",
]
