"""Source-design-only Gaussian GLS oracle contrast, not an empirical rate test.

No function accepts an observed response vector; no random worlds are drawn.
The exact normal distribution assumes BOTH a correct linear conditional mean
and a correctly specified known covariance. Real PROCEED assumptions are NOT
verified by this mathematical benchmark.
"""
from __future__ import annotations
import numpy as np
from scipy.linalg import solve_triangular
from scipy.optimize import brentq
from scipy.stats import norm

ALPHA = 1.0 / 60.0


def cluster_covariance(species, studies):
    """The previous hypothetical DGP covariance, not a fitted real covariance."""
    s=np.asarray(species); t=np.asarray(studies)
    if s.ndim!=1 or t.shape!=s.shape or len(s)==0:
        raise ValueError('Aligned nonempty species and study identifiers required')
    return .25*((s[:,None]==s).astype(float)+(t[:,None]==t).astype(float)+np.eye(len(s)))


def history_contrast(nuisance, history, covariance):
    """Return design-only efficient contrast for a common history coefficient.

    `nuisance` must contain the intended intercept. Collinear controls are
    allowed: only their column span matters. Return weights for verification
    in private memory; the source runner does not export them.
    """
    Z=np.asarray(nuisance,dtype=float)
    h=np.asarray(history,dtype=float)
    S=np.asarray(covariance,dtype=float)
    if (Z.ndim!=2 or h.ndim!=1 or len(h)<2 or Z.shape[0]!=len(h)
        or Z.shape[1]<1 or S.shape!=(len(h),len(h))):
        raise ValueError('Aligned nuisance matrix, history and square covariance required')
    if not (np.isfinite(Z).all() and np.isfinite(h).all() and np.isfinite(S).all()):
        raise ValueError('Finite source geometry and known covariance required')
    if not np.allclose(S,S.T,rtol=1e-12,atol=1e-13):
        raise ValueError('Covariance must be symmetric')
    try:
        L=np.linalg.cholesky(S)
    except np.linalg.LinAlgError as exc:
        raise ValueError('Covariance must be positive definite') from exc
    D=solve_triangular(L,Z,lower=True,check_finite=False)
    g=solve_triangular(L,h,lower=True,check_finite=False)
    scales=np.linalg.norm(D,axis=0)
    nonzero=scales>0
    if not nonzero.any():
        raise ValueError('Nuisance column span is empty; explicit intercept required')
    # Column scaling is an exact span reparameterization, not a ridge penalty.
    normalized=D[:,nonzero]/scales[nonzero]
    U,singular,_=np.linalg.svd(normalized,full_matrices=False)
    rank=int(np.sum(singular>singular[0]*1e-10))
    basis=U[:,:rank]
    r=g-basis@(basis.T@g)
    information=float(r@r)
    if information<=1e-12*max(1.0,float(g@g)):
        raise ValueError('Historical coefficient not identifiable beyond nuisance span')
    weights=solve_triangular(L.T,r/np.sqrt(information),lower=False,check_finite=False)
    annihilation=weights@Z
    scaled_leak=float(np.max(np.abs(annihilation)/np.maximum(1.0,np.linalg.norm(Z,axis=0))))
    variance=float(weights@S@weights)
    if scaled_leak>1e-9 or not np.isclose(variance,1.0,rtol=0,atol=1e-9):
        raise ArithmeticError('Nuisance orthogonality or unit-variance numerical check failed')
    return {
        'information':information,
        'coefficient_standard_error_if_covariance_known':float(1/np.sqrt(information)),
        'nuisance_rank':rank,
        'nuisance_columns':int(Z.shape[1]),
        'scaled_nuisance_leakage':scaled_leak,
        'oracle_null_variance':variance,
        'history_mean_shift_per_unit_coefficient':float(weights@h),
        'z_weights':weights,
    }


def two_sided_power(mean_shift, variance, alpha=ALPHA):
    """Exact two-sided Gaussian rejection probability, not Monte Carlo coverage."""
    if not (np.isfinite(mean_shift) and np.isfinite(variance) and variance>0
            and np.isfinite(alpha) and 0<alpha<1):
        raise ValueError('Finite mean, positive variance and probability alpha required')
    critical=float(norm.isf(alpha/2))
    sd=np.sqrt(variance)
    return float(norm.cdf((-critical-mean_shift)/sd)+norm.sf((critical-mean_shift)/sd))


def effect_for_power(information, target=.8, alpha=ALPHA):
    """Common coefficient needed for target power in the stated oracle model."""
    if not (np.isfinite(information) and information>0 and 0<alpha<target<1):
        raise ValueError('Positive information and alpha < target < 1 required')
    root=brentq(lambda delta:two_sided_power(delta,1.,alpha)-target,0.,20.,xtol=1e-12)
    return float(root/np.sqrt(information))
