"""Post-terminal *oracle design* benchmark; no biological-response input.

Reconstruct original source-only covariates and retain ALL byte checks. A new
mathematical contrast is evaluated analytically, not fitted to rates or
validated by additional random worlds. It cannot authorise old Haldane v1.
"""
from __future__ import annotations
import argparse
import hashlib
import json
from pathlib import Path
import tempfile
import numpy as np
from history_score_oracle import (ALPHA, cluster_covariance, history_contrast,
                                   two_sided_power, effect_for_power)

HERE=Path(__file__).resolve().parent
CLOSED='TERMINAL_NOT_EVALUABLE_H_SYNTHETIC_NULL_CEILING_FAIL_NO_BIOLOGICAL_OUTCOME_OPENING'


def digest(a):
    return hashlib.sha256(np.asarray(a,dtype='<f8').tobytes(order='C')).hexdigest()


def analyse_design(X,H,species,studies,roles):
    X=np.asarray(X,float); H=np.asarray(H,float)
    s=np.asarray(species); st=np.asarray(studies); role=np.asarray(roles)
    if X.shape!=(165,21) or H.shape!=(165,) or s.shape!=(165,) or st.shape!=(165,) or role.shape!=(165,):
        raise ValueError('Original source predictor dimensions changed')
    if not np.isfinite(X).all() or not np.isfinite(H).all():
        raise ValueError('Finite source covariates required')
    train=role=='train'; held=role=='heldout'
    if (int(train.sum()),int(held.sum()))!=(91,74):
        raise ValueError('Original pair-role assignment changed')
    if (len(np.unique(s[train])),len(np.unique(s[held])))!=(29,27):
        raise ValueError('Original source species-role counts changed')
    if set(s[train])&set(s[held]) or set(st[train])&set(st[held]):
        raise ValueError('Source train/heldout study or species leakage')
    w=np.zeros(165)
    for sp in np.unique(s[train]):
        mask=train&(s==sp);w[mask]=1/(29*int(mask.sum()))
    hmean=float(w@H); hsd=float(np.sqrt(w@((H-hmean)**2)))
    if hsd<=1e-12:raise ValueError('No original training history variance')
    h=(H[held]-hmean)/hsd
    # Scaling changes no nuisance span. This is not coefficient shrinkage.
    xmean=w@X
    xsd=np.sqrt(w@((X-xmean)**2));xsd=np.where(xsd>1e-10,xsd,1.)
    Z=np.column_stack([np.ones(int(held.sum())),(X[held]-xmean)/xsd])
    true_sigma=cluster_covariance(s[held],st[held])
    full=history_contrast(Z,h,true_sigma)
    u=full['z_weights']
    variance=float(u@true_sigma@u)
    shift=float(u@h)
    diagonal=history_contrast(Z,h,np.diag(np.diag(true_sigma)))
    ud=diagonal['z_weights'];actual_diag_variance=float(ud@true_sigma@ud)
    common={k:v for k,v in full.items() if k!='z_weights'}
    common.update({
        'analytic_null_rejection_probability':two_sided_power(0,variance,ALPHA),
        'power_grid':[{'gamma_per_train_history_SD':g,
                       'analytic_power_under_known_gaussian_covariance':two_sided_power(g*shift,variance,ALPHA)}
                      for g in (0.,.25,.5,1.)],
        'gamma_for_80pct_oracle_power':effect_for_power(full['information'],.8,ALPHA),
    })
    return {
        'heldout_rows':int(held.sum()),'heldout_species':27,
        'heldout_original_studies':len(np.unique(st[held])),
        'training_history_equal_species_SD':hsd,
        'training_history_equal_species_mean':hmean,
        'two_sided_alpha':ALPHA,
        'known_covariance_oracle':common,
        'prespecified_diagonal_assumption_sensitivity':{
            'assumed_null_variance':diagonal['oracle_null_variance'],
            'actual_variance_under_same_clustered_gaussian_model':actual_diag_variance,
            'analytic_null_rejection_probability_under_same_clustered_model':two_sided_power(0,actual_diag_variance,ALPHA),
            'not_used_to_choose_a_winning_method':True,
        },
        'empirical_response_opening_authorized':False,
        'oracle_is_real_biological_validation':False,
        'needs_correct_known_covariance_and_correct_mean_model':True,
        'common_coefficient_assumed_not_established':True,
        'not_out_of_species_predictive_gain':True,
    }


def main():
    p=argparse.ArgumentParser();p.add_argument('--output',required=True,type=Path)
    args=p.parse_args()
    receipt_path=HERE/'proceed70_haldane_synthetic_terminal_receipt_v01.json'
    original_bytes=receipt_path.read_bytes()
    receipt=json.loads(original_bytes)
    if receipt['status']!=CLOSED:raise ValueError('Old terminal closure must remain in force')
    output={
        'schema':'ttf_history_conditional_signal_oracle_design_benchmark_v01',
        'old_terminal_status_unchanged':CLOSED,
        'empirical_rate_or_fitness_values_read':False,
        'additional_null_or_positive_worlds_drawn':0,
        'old_thresholds_model_or_hashes_changed':False,
        'no_empirical_opening_authorized':True,
        'method_assumptions':'Hypothetical Gaussian linear common-effect model with EXACT KNOWN covariance; not estimated from biological data',
    }
    args.output.parent.mkdir(parents=True,exist_ok=True)
    try:
        from terminal_numeric_provenance import code_identity,prepare_sources,worker
        output['source_code_identity']=code_identity()
        if not all(v['identical'] for v in output['source_code_identity'].values()):
            raise ValueError('Source helpers changed since closed original calibration')
        with tempfile.TemporaryDirectory(prefix='ttf-private-oracle-') as td:
            work=Path(td)
            output['source_identity']=prepare_sources(work,receipt)
            worker(work,'oracle_source')
            lock=receipt['source_pre_registered']
            with np.load(work/'oracle_source.npz',allow_pickle=False) as a:
                xsha=digest(a['X']);hsha=digest(a['H'])
                output['source_X_sha256']=xsha;output['source_H_sha256']=hsha
                if (xsha!=lock['source_21_baseline_numeric_SHA256'] or
                    hsha!=lock['source_one_history_numeric_SHA256']):
                    raise ValueError('Original X/H bytes do not match; oracle design not evaluated')
                output['oracle_design']=analyse_design(a['X'],a['H'],a['sp'],a['st'],a['role'])
        if receipt_path.read_bytes()!=original_bytes:
            raise ValueError('Original terminal receipt changed')
        output['status']='ORACLE_DESIGN_BENCHMARK_COMPLETE_OLD_BIOLOGICAL_STOP_UNCHANGED'
    except Exception as exc:
        output['status']='ORACLE_DESIGN_DIAGNOSTIC_STOP_OLD_BIOLOGICAL_STOP_UNCHANGED'
        output['error']=type(exc).__name__+': '+str(exc)[:500]
        args.output.write_text(json.dumps(output,sort_keys=True,indent=2)+'\n')
        raise
    args.output.write_text(json.dumps(output,sort_keys=True,indent=2)+'\n')
    print(json.dumps(output,sort_keys=True))

if __name__=='__main__':main()
