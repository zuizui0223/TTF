"""Post-terminal unknown-variance DESIGN benchmark, never original outcomes.

All source byte checks remain mandatory. Species/study intercepts are removed
as nuisance directions, not assigned known variances. This is not an empirical
reopening or a replacement calibration for the failed Haldane test.
"""
from __future__ import annotations
import argparse
import hashlib
import json
from pathlib import Path
import tempfile
import numpy as np
from scipy.stats import t
from history_absorbed_design import absorbed_design,student_power,detectable_effect,ALPHA

HERE=Path(__file__).resolve().parent
CLOSED='TERMINAL_NOT_EVALUABLE_H_SYNTHETIC_NULL_CEILING_FAIL_NO_BIOLOGICAL_OUTCOME_OPENING'


def digest(a):
    return hashlib.sha256(np.asarray(a,dtype='<f8').tobytes(order='C')).hexdigest()


def matches_original(X,H,expected_X,expected_H):
    return digest(X)==expected_X and digest(H)==expected_H


def analyse_source(X,H,species,studies,roles):
    X=np.asarray(X,float);H=np.asarray(H,float)
    sp=np.asarray(species);st=np.asarray(studies);role=np.asarray(roles)
    if (X.shape!=(165,21) or H.shape!=(165,) or sp.shape!=(165,) or
        st.shape!=(165,) or role.shape!=(165,) or not np.isfinite(X).all() or not np.isfinite(H).all()):
        raise ValueError('Original finite source X/H and labels changed')
    train=role=='train';held=role=='heldout'
    if (int(train.sum()),int(held.sum()))!=(91,74):raise ValueError('Original roles changed')
    if (len(np.unique(sp[train])),len(np.unique(sp[held])))!=(29,27):raise ValueError('Original species counts changed')
    if (len(np.unique(st[train])),len(np.unique(st[held])))!=(48,45):raise ValueError('Original study counts changed')
    if set(sp[train])&set(sp[held]) or set(st[train])&set(st[held]):raise ValueError('Source split leaks studies/species')
    w=np.zeros(165)
    for s in np.unique(sp[train]):
        at=train&(sp==s);w[at]=1/(29*int(at.sum()))
    hmean=float(w@H);hsd=float(np.sqrt(w@((H-hmean)**2)))
    if hsd<=1e-12:raise ValueError('No source training history variation')
    xmean=w@X;xsd=np.sqrt(w@((X-xmean)**2));xsd=np.where(xsd>1e-10,xsd,1.)
    Z=np.column_stack([np.ones(74),(X[held]-xmean)/xsd])
    h=(H[held]-hmean)/hsd
    d=absorbed_design(Z,h,sp[held],st[held])
    public={k:v for k,v in d.items() if not isinstance(v,np.ndarray)}
    if d['status']=='IDENTIFIABLE_UNKNOWN_SCALE_GAUSSIAN_DESIGN':
        df=d['residual_df'];energy=d['history_residual_energy'];u=d['unit_contrast']
        shares=np.array([np.sum(u[sp[held]==s]**2) for s in np.unique(sp[held])])
        public.update({
          'conditional_gaussian_null_rejection_probability':student_power(0,energy,df,.5),
          'student_t_two_sided_cutoff':float(t.isf(ALPHA/2,df)),
          'hypothetical_residual_sigma_for_power_only':.5,
          'coefficient_SE_multiplier_of_unknown_sigma':float(1/np.sqrt(energy)),
          'power_grid':[{'gamma_per_training_history_SD':g,
                         'power_if_iid_residual_sigma_point5':student_power(g,energy,df,.5)} for g in (0.,.25,.5,1.)],
          'gamma_80pct_power_if_residual_sigma_point5':detectable_effect(energy,df,.5),
          'effective_species_contrast_energy':float(1/np.sum(shares**2)),
          'maximum_species_contrast_energy_share':float(shares.max()),
          'nuisance_group_variance_values_required':False,
          'residual_variance_known_required':False,
        })
    return {'heldout_rows':74,'heldout_species':27,'heldout_studies':45,
            'training_history_SD':hsd,'two_sided_alpha':ALPHA,
            'within_group_design':public,
            'requires_iid_homoscedastic_normal_residuals_after_group_offsets':True,
            'handles_arbitrary_unknown_additive_species_study_offsets':True,
            'does_not_handle_unknown_random_slopes_or_spatial_phylogenetic_residuals':True,
            'removes_between_group_history_and_changes_estimand':True,
            'common_history_coefficient_assumed_not_established':True,
            'real_phenotype_validity_established':False,
            'empirical_opening_authorized':False,
            'not_unknown_species_prediction_gain_or_causal_adaptation':True}


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--output',required=True,type=Path)
    args=ap.parse_args()
    path=HERE/'proceed70_haldane_synthetic_terminal_receipt_v01.json'
    before=path.read_bytes();terminal=json.loads(before)
    if terminal['status']!=CLOSED:raise ValueError('Old terminal closure required')
    out={'schema':'ttf_history_absorbed_unknown_variance_design_v01',
         'old_terminal_status_unchanged':CLOSED,
         'actual_phenotype_or_fitness_values_read':False,
         'new_random_worlds_drawn':0,'old_sources_parameters_or_thresholds_changed':False,
         'empirical_response_opening_authorized':False}
    args.output.parent.mkdir(parents=True,exist_ok=True)
    try:
        from terminal_numeric_provenance import code_identity,prepare_sources,worker
        out['source_code_identity']=code_identity()
        if not all(v['identical'] for v in out['source_code_identity'].values()):
            raise ValueError('Frozen source-processing code modified')
        with tempfile.TemporaryDirectory(prefix='ttf-private-absorbed-') as td:
            work=Path(td);out['source_identity']=prepare_sources(work,terminal)
            worker(work,'unchanged_source')
            lock=terminal['source_pre_registered']
            with np.load(work/'unchanged_source.npz',allow_pickle=False) as data:
                out['source_X_sha256']=digest(data['X']);out['source_H_sha256']=digest(data['H'])
                if not matches_original(data['X'],data['H'],lock['source_21_baseline_numeric_SHA256'],lock['source_one_history_numeric_SHA256']):
                    raise ValueError('Original source X/H failed exact-byte identity')
                out['design']=analyse_source(data['X'],data['H'],data['sp'],data['st'],data['role'])
        if path.read_bytes()!=before:raise ValueError('Original terminal decision was modified')
        out['status']='UNKNOWN_VARIANCE_DESIGN_DIAGNOSTIC_COMPLETE_OLD_STOP_UNCHANGED'
    except Exception as exc:
        out.update(status='DESIGN_DIAGNOSTIC_STOP_OLD_STOP_UNCHANGED',error=type(exc).__name__+': '+str(exc)[:450])
        args.output.write_text(json.dumps(out,indent=2,sort_keys=True)+'\n')
        raise
    args.output.write_text(json.dumps(out,indent=2,sort_keys=True)+'\n')
    print(json.dumps(out,sort_keys=True))

if __name__=='__main__':main()
