from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .relational_environment import canonical_pca, project_whiten, schoener_d


@dataclass(frozen=True)
class CommonClimateSpace:
    mean: np.ndarray
    sd: np.ndarray
    eigval: np.ndarray
    eigvec: np.ndarray
    axes: int


def fit_common_climate_space(
    pooled_current_environment: np.ndarray,
    *,
    axes: int = 4,
) -> CommonClimateSpace:
    x=np.asarray(pooled_current_environment,dtype=float)
    if x.ndim!=2 or x.shape[1]<axes:
        raise ValueError("pooled climate matrix has fewer columns than requested axes")
    mean,sd,eigval,eigvec=canonical_pca(x)
    if axes<1 or axes>eigvec.shape[1]:
        raise ValueError("invalid retained climate-space axes")
    return CommonClimateSpace(
        mean=np.asarray(mean,float),
        sd=np.asarray(sd,float),
        eigval=np.asarray(eigval,float),
        eigvec=np.asarray(eigvec,float),
        axes=int(axes),
    )


def whiten_environment(
    environment: np.ndarray,
    space: CommonClimateSpace,
) -> np.ndarray:
    return project_whiten(
        environment,
        space.mean,
        space.sd,
        space.eigval,
        space.eigvec,
        axes=space.axes,
    )


def scott_factor(n: int, dimensions: int) -> float:
    if n<2 or dimensions<1:
        raise ValueError("Scott bandwidth requires n>=2 and dimensions>=1")
    return float(n ** (-1.0 / (dimensions + 4.0)))


def gaussian_kde_density(
    training_points: np.ndarray,
    query_points: np.ndarray,
    *,
    chunk_size: int = 4096,
) -> np.ndarray:
    x=np.asarray(training_points,dtype=float)
    q=np.asarray(query_points,dtype=float)
    if x.ndim!=2 or q.ndim!=2 or x.shape[1]!=q.shape[1]:
        raise ValueError("training/query climate points must be n x d with equal d")
    if len(x)<2 or len(q)<1 or not np.isfinite(x).all() or not np.isfinite(q).all():
        raise ValueError("KDE climate points must be finite and non-empty")
    if chunk_size<1:
        raise ValueError("chunk_size must be positive")

    h=scott_factor(len(x),x.shape[1])
    out=np.empty(len(q),dtype=float)
    inv=1.0/(h*h)
    for start in range(0,len(q),int(chunk_size)):
        stop=min(len(q),start+int(chunk_size))
        delta=q[start:stop,None,:]-x[None,:,:]
        sq=np.sum(delta*delta,axis=2)
        out[start:stop]=np.exp(-0.5*sq*inv).mean(axis=1)
    return out


def normalized_suitability(
    training_current_environment: np.ndarray,
    grid_environment: np.ndarray,
    space: CommonClimateSpace,
    *,
    chunk_size: int = 4096,
) -> np.ndarray:
    training=whiten_environment(training_current_environment,space)
    grid=whiten_environment(grid_environment,space)
    density=gaussian_kde_density(training,grid,chunk_size=chunk_size)
    total=float(density.sum())
    if not np.isfinite(total) or total<=0:
        raise ValueError("climatic-envelope suitability has zero/non-finite mass")
    return density/total


def extrapolation_fraction(
    grid_environment: np.ndarray,
    pooled_current_environment: np.ndarray,
    space: CommonClimateSpace,
) -> float:
    grid=whiten_environment(grid_environment,space)
    pooled=whiten_environment(pooled_current_environment,space)
    low=pooled.min(axis=0)
    high=pooled.max(axis=0)
    outside=np.any((grid<low[None,:]) | (grid>high[None,:]),axis=1)
    return float(np.mean(outside))


def current_lgm_suitability(
    species_current_environment: np.ndarray,
    current_grid_environment: np.ndarray,
    lgm_grid_environment: np.ndarray,
    space: CommonClimateSpace,
    *,
    chunk_size: int = 4096,
) -> tuple[np.ndarray,np.ndarray]:
    current=normalized_suitability(
        species_current_environment,
        current_grid_environment,
        space,
        chunk_size=chunk_size,
    )
    lgm=normalized_suitability(
        species_current_environment,
        lgm_grid_environment,
        space,
        chunk_size=chunk_size,
    )
    return current,lgm


def climatic_relation(
    left_suitability: np.ndarray,
    right_suitability: np.ndarray,
) -> float:
    return schoener_d(left_suitability,right_suitability)
