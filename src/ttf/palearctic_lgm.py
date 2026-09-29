from __future__ import annotations

from dataclasses import dataclass
import hashlib

import numpy as np


ROLE_NAMESPACE = "palearctic-insect-lgm-role-v0.1"
ARCHIVE_SHA256 = "5a0fd9ac25893c749d14186fbcce4a46b99163c9d810b36e40eebce61a5"


def palearctic_role_order(species: list[str]) -> list[str]:
    names=list(map(str,species))
    if len(names) != len(set(names)):
        raise ValueError("species names must be unique")
    return sorted(
        names,
        key=lambda name: (
            hashlib.sha256(
                f"{ROLE_NAMESPACE}|{ARCHIVE_SHA256}|{name}".encode("utf-8")
            ).hexdigest(),
            name,
        ),
    )


@dataclass(frozen=True)
class ClimateScaling:
    mean: np.ndarray
    sd: np.ndarray
    lower_z: np.ndarray
    upper_z: np.ndarray


@dataclass(frozen=True)
class RidgeLogisticFit:
    coefficients: np.ndarray
    converged: bool
    iterations: int
    objective: float


def target_group_background_cells(
    lon_lat: np.ndarray,
    *,
    resolution_degrees: float = 0.25,
    maximum_cells: int = 20_000,
    namespace: str = "palearctic-insect-lgm-background-v0.2",
) -> np.ndarray:
    """Return deterministic unique target-group background cell centres.

    Input coordinates must already have passed the formal Palearctic realm mask.
    The function has no access to any genetic response.
    """
    x=np.asarray(lon_lat,dtype=float)
    if x.ndim != 2 or x.shape[1] != 2 or not len(x) or not np.isfinite(x).all():
        raise ValueError("lon_lat must be a non-empty finite n x 2 array")
    if resolution_degrees <= 0 or maximum_cells < 1:
        raise ValueError("invalid background grid settings")

    lon=x[:,0]
    lat=x[:,1]
    if np.any((lon < -180) | (lon > 180) | (lat < -90) | (lat > 90)):
        raise ValueError("invalid longitude/latitude")

    # Canonical global cell origin is (-180,-90). Use integer cell indices to
    # avoid float-string ambiguity in deduplication and ranking.
    ix=np.floor((lon + 180.0) / resolution_degrees).astype(np.int64)
    iy=np.floor((lat + 90.0) / resolution_degrees).astype(np.int64)
    nx=int(round(360.0 / resolution_degrees))
    ny=int(round(180.0 / resolution_degrees))
    ix=np.clip(ix,0,nx-1)
    iy=np.clip(iy,0,ny-1)
    pairs=sorted(set(zip(ix.tolist(),iy.tolist())))

    def key(pair: tuple[int,int]) -> tuple[str,int,int]:
        i,j=pair
        centre_lon=-180.0+(i+0.5)*resolution_degrees
        centre_lat=-90.0+(j+0.5)*resolution_degrees
        token=f"{namespace}|{centre_lon:.6f}|{centre_lat:.6f}"
        return hashlib.sha256(token.encode("utf-8")).hexdigest(),i,j

    selected=sorted(pairs,key=key)[:maximum_cells]
    return np.asarray([
        [
            -180.0+(i+0.5)*resolution_degrees,
            -90.0+(j+0.5)*resolution_degrees,
        ]
        for i,j in selected
    ],dtype=float)


def fit_climate_scaling(background_current: np.ndarray) -> ClimateScaling:
    x=np.asarray(background_current,dtype=float)
    if x.ndim != 2 or x.shape[0] < 20 or x.shape[1] != 4 or not np.isfinite(x).all():
        raise ValueError("background_current must be finite n x 4 with n>=20")
    mean=x.mean(axis=0)
    sd=x.std(axis=0,ddof=0)
    if np.any(sd <= np.finfo(float).eps):
        raise ValueError("background climate has a constant variable")
    z=(x-mean)/sd
    lower=np.percentile(z,0.5,axis=0)
    upper=np.percentile(z,99.5,axis=0)
    if np.any(upper <= lower):
        raise ValueError("degenerate climate clamp")
    return ClimateScaling(mean=mean,sd=sd,lower_z=lower,upper_z=upper)


def climate_features(values: np.ndarray, scaling: ClimateScaling) -> np.ndarray:
    x=np.asarray(values,dtype=float)
    if x.ndim != 2 or x.shape[1] != 4 or not np.isfinite(x).all():
        raise ValueError("climate values must be finite n x 4")
    z=(x-scaling.mean)/scaling.sd
    z=np.clip(z,scaling.lower_z,scaling.upper_z)
    return np.column_stack((z,z*z))


def _sigmoid(x: np.ndarray) -> np.ndarray:
    out=np.empty_like(x,dtype=float)
    positive=x>=0
    out[positive]=1.0/(1.0+np.exp(-x[positive]))
    expx=np.exp(x[~positive])
    out[~positive]=expx/(1.0+expx)
    return out


def _objective(
    X: np.ndarray,
    y: np.ndarray,
    weights: np.ndarray,
    beta: np.ndarray,
    ridge_lambda: float,
) -> float:
    eta=X@beta
    # log(1+exp(eta))-y*eta, stable.
    loss=np.logaddexp(0.0,eta)-y*eta
    penalty=0.5*ridge_lambda*float(np.dot(beta[1:],beta[1:]))
    return float(np.dot(weights,loss)+penalty)


def fit_weighted_ridge_logistic(
    presence_features: np.ndarray,
    background_features: np.ndarray,
    *,
    ridge_lambda: float = 1.0,
    presence_total_weight: float = 0.5,
    background_total_weight: float = 0.5,
    max_iter: int = 2000,
    tolerance: float = 1e-10,
) -> RidgeLogisticFit:
    """Fit deterministic L2-regularized presence-background logistic model.

    The intercept is not penalized. Presence and background class weights are
    normalized separately to the frozen total weights.
    """
    p=np.asarray(presence_features,dtype=float)
    b=np.asarray(background_features,dtype=float)
    if (
        p.ndim != 2 or b.ndim != 2 or p.shape[1] != b.shape[1]
        or len(p) < 2 or len(b) < 2 or not np.isfinite(p).all()
        or not np.isfinite(b).all()
    ):
        raise ValueError("presence/background features must be finite compatible matrices")
    if ridge_lambda <= 0 or presence_total_weight <= 0 or background_total_weight <= 0:
        raise ValueError("invalid logistic regularization/weights")

    X=np.vstack((p,b))
    X=np.column_stack((np.ones(len(X)),X))
    y=np.r_[np.ones(len(p)),np.zeros(len(b))]
    weights=np.r_[
        np.full(len(p),presence_total_weight/len(p)),
        np.full(len(b),background_total_weight/len(b)),
    ]

    beta=np.zeros(X.shape[1],dtype=float)
    # Balanced frozen weights imply zero initial intercept.
    penalty=np.eye(X.shape[1],dtype=float)*ridge_lambda
    penalty[0,0]=0.0
    last=_objective(X,y,weights,beta,ridge_lambda)
    converged=False

    for iteration in range(1,max_iter+1):
        mu=_sigmoid(X@beta)
        grad=X.T@(weights*(mu-y))+penalty@beta
        curvature=weights*mu*(1.0-mu)
        H=X.T@(X*curvature[:,None])+penalty
        try:
            step=np.linalg.solve(H,grad)
        except np.linalg.LinAlgError:
            step=np.linalg.lstsq(H,grad,rcond=None)[0]

        # Deterministic backtracking preserves objective descent.
        scale=1.0
        accepted=False
        for _ in range(40):
            trial=beta-scale*step
            value=_objective(X,y,weights,trial,ridge_lambda)
            if value <= last + 1e-14:
                beta=trial
                accepted=True
                break
            scale*=0.5
        if not accepted:
            break

        delta=float(np.max(np.abs(scale*step)))
        last=value
        if delta <= tolerance:
            converged=True
            break

    return RidgeLogisticFit(
        coefficients=beta,
        converged=converged,
        iterations=iteration,
        objective=last,
    )


def predict_relative_suitability(features: np.ndarray, fit: RidgeLogisticFit) -> np.ndarray:
    x=np.asarray(features,dtype=float)
    if x.ndim != 2 or x.shape[1]+1 != len(fit.coefficients) or not np.isfinite(x).all():
        raise ValueError("prediction features do not match fitted model")
    eta=fit.coefficients[0]+x@fit.coefficients[1:]
    # Odds/rate scale is the natural relative-intensity scale for a
    # presence-background logistic fit; normalize to unit mass afterwards.
    eta=eta-float(np.max(eta))
    score=np.exp(eta)
    total=float(score.sum())
    if not np.isfinite(total) or total <= 0:
        raise ValueError("non-finite suitability")
    return score/total


def schoener_d(left: np.ndarray, right: np.ndarray) -> float:
    a=np.asarray(left,dtype=float)
    b=np.asarray(right,dtype=float)
    if a.ndim != 1 or b.shape != a.shape or not np.isfinite(a).all() or not np.isfinite(b).all():
        raise ValueError("suitability vectors must be finite equal-length vectors")
    if not np.isclose(a.sum(),1.0,atol=1e-10) or not np.isclose(b.sum(),1.0,atol=1e-10):
        raise ValueError("suitability vectors must sum to one")
    return float(np.clip(1.0-0.5*np.abs(a-b).sum(),0.0,1.0))


def relation_matrix(suitability: np.ndarray) -> np.ndarray:
    s=np.asarray(suitability,dtype=float)
    if s.ndim != 2 or s.shape[0] < 2 or not np.isfinite(s).all():
        raise ValueError("suitability must be species x cells")
    if not np.allclose(s.sum(axis=1),1.0,atol=1e-10):
        raise ValueError("each suitability surface must sum to one")
    n=len(s)
    out=np.eye(n,dtype=float)
    for i in range(n):
        for j in range(i+1,n):
            value=schoener_d(s[i],s[j])
            out[i,j]=value
            out[j,i]=value
    return out
