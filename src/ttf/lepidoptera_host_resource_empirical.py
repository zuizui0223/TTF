from __future__ import annotations

from dataclasses import dataclass
from typing import Mapping

import numpy as np

from .lepidoptera_host_resource_qualification_v02 import nuisance_envelope_pvalues
from .lepidoptera_trait_gradient_v02 import target_equal_mean_correlation


@dataclass(frozen=True)
class EmpiricalHostResourceScore:
    statistic: float
    p_value: float
    positive: bool
    per_target: Mapping[int, float]


def frozen_pair_surface(design_npz):
    names=tuple(map(str,design_npz["species_order"]))
    target=np.asarray(design_npz["target_index"],np.int64)
    source=np.asarray(design_npz["source_index"],np.int64)
    residual=np.asarray(design_npz["host_resource_residual"],float)
    lengths=np.asarray(design_npz["alignment_length"],np.int64)
    offsets=np.asarray(design_npz["alignment_offset"],np.int64)
    expected=np.r_[0,np.cumsum(lengths)[:-1]]
    if not np.array_equal(offsets,expected):
        raise RuntimeError("host-resource alignment offset drift")
    pair_id=np.repeat(np.arange(len(target),dtype=np.int64),lengths)
    return names,target,source,residual,pair_id


def pair_transfer_congruence(
    response_by_species: Mapping[str,np.ndarray],
    species_order: tuple[str,...],
    target_index: np.ndarray,
    source_index: np.ndarray,
    pair_id: np.ndarray,
    alignment_target: np.ndarray,
    alignment_source: np.ndarray,
) -> np.ndarray:
    if set(response_by_species)!=set(species_order):
        raise ValueError("response species differ from frozen host-resource design")
    edge_counts=np.asarray(
        [len(np.asarray(response_by_species[name],float)) for name in species_order],
        np.int64,
    )
    edge_offsets=np.r_[0,np.cumsum(edge_counts)]
    flat=np.concatenate([
        np.asarray(response_by_species[name],float) for name in species_order
    ])
    gt=edge_offsets[target_index][pair_id]+np.asarray(alignment_target,np.int64)
    gs=edge_offsets[source_index][pair_id]+np.asarray(alignment_source,np.int64)
    x=flat[gt]; y=flat[gs]
    n_pairs=len(target_index)
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
    return out


def score_empirical_host_resource(
    response_by_species: Mapping[str,np.ndarray],
    design_npz,
    *,
    private_reference: np.ndarray,
    geometry_reference: np.ndarray,
    breadth_reference: np.ndarray,
    alpha: float=.05,
) -> EmpiricalHostResourceScore:
    names,target,source,residual,pair_id=frozen_pair_surface(design_npz)
    pair_score=pair_transfer_congruence(
        response_by_species,names,target,source,pair_id,
        np.asarray(design_npz["alignment_target"],np.int64),
        np.asarray(design_npz["alignment_source"],np.int64),
    )
    statistic,per_target=target_equal_mean_correlation(pair_score,residual,target)
    p=float(nuisance_envelope_pvalues(
        np.asarray([statistic],float),
        private_reference=np.asarray(private_reference,float),
        geometry_reference=np.asarray(geometry_reference,float),
        breadth_reference=np.asarray(breadth_reference,float),
    )[0])
    return EmpiricalHostResourceScore(
        statistic=float(statistic),
        p_value=p,
        positive=bool(float(statistic)>0 and p<=float(alpha)),
        per_target=per_target,
    )


__all__=[
    "EmpiricalHostResourceScore",
    "frozen_pair_surface",
    "pair_transfer_congruence",
    "score_empirical_host_resource",
]
