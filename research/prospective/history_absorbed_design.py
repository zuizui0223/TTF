"""Source-only Gaussian design with unknown species/study offsets and iid scale.

No biological outcome argument and no random worlds. Group intercepts are
absorbed rather than assigned known covariance components. Student-t validity
still requires a correct common linear effect and iid normal residuals.
"""
from __future__ import annotations
import numpy as np
from scipy.optimize import brentq
from scipy.stats import nct,t

ALPHA=1/60
RANK_RTOL=1e-10


def incidence(labels):
    a=np.asarray(labels)
    if a.ndim!=1 or len(a)<3 or any(v is None or not str(v).strip() for v in a):
        raise ValueError('Nonempty aligned original group labels required')
    a=a.astype(str)
    return (a[:,None]==np.unique(a)[None,:]).astype(float)


def basis(matrix):
    """Fixed column-normalized rank rule; invariant to benign re-scalings."""
    x=np.asarray(matrix,float)
    if x.ndim!=2 or not np.isfinite(x).all():
        raise ValueError('Finite two-dimensional nuisance matrix required')
    norms=np.linalg.norm(x,axis=0)
    use=norms>0
    if not np.any(use):return np.zeros((len(x),0))
    u,s,_=np.linalg.svd(x[:,use]/norms[use],full_matrices=False)
    return u[:,s>RANK_RTOL*s[0]]


def absorbed_design(Z,h,species,studies):
    """Return design operators only. Never accepts phenotype y or source rates."""
    Z=np.asarray(Z,float);h=np.asarray(h,float)
    if (Z.ndim!=2 or h.shape!=(len(Z),) or len(Z)<3 or
        not np.isfinite(Z).all() or not np.isfinite(h).all()):
        raise ValueError('Finite aligned source controls and historical predictor required')
    A=incidence(species);B=incidence(studies)
    if len(A)!=len(Z) or len(B)!=len(Z):raise ValueError('Group rows misaligned')
    base=np.column_stack([np.ones(len(Z)),Z])
    W=np.column_stack([base,A,B])
    Q=basis(W);Qbase=basis(base);Qgroup=basis(np.column_stack([A,B]))
    r=h-Q@(Q.T@h)
    rb=h-Qbase@(Qbase.T@h)
    energy=float(r@r);base_energy=float(rb@rb)
    df=len(h)-Q.shape[1]-1
    out={'n':len(h),'nuisance_rank':Q.shape[1],
         'base_nuisance_rank':Qbase.shape[1],'group_incidence_rank':Qgroup.shape[1],
         'residual_df':df,'history_residual_energy':energy,
         'base_only_history_energy':base_energy,
         'history_energy_retention_after_group_absorption':energy/base_energy if base_energy>1e-20 else None,
         'rank_relative_tolerance':RANK_RTOL}
    if energy<=RANK_RTOL**2*max(1.,float(h@h)):
        out['status']='NO_WITHIN_GROUP_HISTORY_INFORMATION'
        return out
    if df<=0:
        out['status']='NO_RESIDUAL_VARIANCE_DF'
        return out
    u=r/np.sqrt(energy)
    M=np.eye(len(h))-Q@Q.T-np.outer(u,u)
    out.update(status='IDENTIFIABLE_UNKNOWN_SCALE_GAUSSIAN_DESIGN',
               unit_contrast=u,residual_projector=(M+M.T)/2,nuisance_matrix=W,
               max_nuisance_leakage=float(np.max(np.abs(u@W))),
               unit_variance=float(u@u),
               projection_numerator_denominator_orthogonality=float(np.linalg.norm(M@u)))
    return out


def student_power(gamma,energy,df,sigma=.5):
    """Conditional Gaussian power with unknown residual scale, not observed power."""
    if (not all(np.isfinite(v) for v in (gamma,energy,df,sigma)) or
        energy<=0 or df<1 or int(df)!=df or sigma<=0):
        raise ValueError('Positive information, integer residual df and noise scale required')
    cutoff=t.isf(ALPHA/2,int(df))
    delta=float(gamma)*np.sqrt(float(energy))/float(sigma)
    return float(nct.cdf(-cutoff,int(df),delta)+nct.sf(cutoff,int(df),delta))


def detectable_effect(energy,df,sigma=.5):
    student_power(0,energy,df,sigma)
    upper=1.
    while student_power(upper,energy,df,sigma)<.8:
        upper*=2
        if upper>1e12:raise ValueError('Detectable effect outside finite design search')
    return float(brentq(lambda g:student_power(g,energy,df,sigma)-.8,0.,upper,xtol=1e-12))
