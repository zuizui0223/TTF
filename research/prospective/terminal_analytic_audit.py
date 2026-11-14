"""Post-terminal analytical diagnosis ONLY: no extra null worlds or real rates.

Reconstructs unchanged source geometry and integrates the already frozen
Gaussian simulation model analytically. No hypothesis is reopened, no model
is tuned, and no empirical outcome is read. Output is aggregate diagnostics.
"""
from __future__ import annotations
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np

ALPHA=.1
FROZEN_X='e4d2f7e5b0ebf900018ae9f8ed87715a90101181b565c83b4215290f6fa07d84'
FROZEN_H='b64289ab8a6719a5655ef3f585887c659cd3441455ffff9f92237b2f548f4c9e'
FROZEN_CLIMATE='7d03b2f158ade223b752a2b7d6ac57c38fd163e0b44f5566a01e8e6fe6d6b4fb'
FROZEN_COORD='aed6486ee9aa49b0d454426c4d5964ebd79a8a96db2aac452c80fd5b7f52af32'


def wilson(k,n):
    if n<1 or not 0<=k<=n: raise ValueError('invalid binomial counts')
    z=1.959963984540054;p=k/n;d=1+z*z/n
    c=(p+z*z/(2*n))/d
    s=z*np.sqrt(p*(1-p)/n+z*z/(4*n*n))/d
    return [float(max(0,c-s)),float(min(1,c+s))]


def scale(X,w):
    mu=w@X;sd=np.sqrt(np.maximum(w@((X-mu)**2),0))
    return mu,np.where(sd>1e-10,sd,1.)


def ridge_operator(X,train,w):
    """Linear map from training outcomes to heldout predictions, frozen alpha."""
    X=np.asarray(X,float);train=np.asarray(train,bool);w=np.asarray(w,float)
    if X.ndim!=2 or len(X)!=len(train) or len(w)!=train.sum():
        raise ValueError('invalid ridge geometry')
    if not np.isfinite(X).all() or np.any(w<=0) or not np.isclose(w.sum(),1):
        raise ValueError('finite predictors and normalized positive weights required')
    mu,sd=scale(X[train],w);Z=(X-mu)/sd
    C=np.eye(train.sum())-np.ones((train.sum(),1))@w[None,:]
    gram=Z[train].T@(w[:,None]*Z[train])+ALPHA*np.eye(Z.shape[1])
    coeff=np.linalg.solve(gram,Z[train].T*w[None,:])@C
    return np.ones(((~train).sum(),1))@w[None,:]+Z[~train]@coeff


def noise_covariance(species,studies):
    sp=np.asarray(species);st=np.asarray(studies)
    if len(sp)!=len(st): raise ValueError('misaligned species and study')
    return .25*((sp[:,None]==sp).astype(float)+(st[:,None]==st).astype(float)+np.eye(len(sp)))


def error_operator(pred,train):
    train=np.asarray(train,bool);idx=np.flatnonzero(~train)
    A=np.zeros((len(idx),len(train)))
    A[:,train]=-pred;A[np.arange(len(idx)),idx]=1
    return A


def square_covariance(mu1,mu2,cross_cov):
    return 2*cross_cov**2+4*np.outer(mu1,mu2)*cross_cov


def loss_difference_moments(A,B,mean,cov):
    """Gaussian exact first/second moments of e_base^2-e_aug^2 by test row."""
    u=A@mean;v=B@mean
    aa=A@cov@A.T;bb=B@cov@B.T;ab=A@cov@B.T
    expected=u*u+np.diag(aa)-v*v-np.diag(bb)
    variance=(square_covariance(u,u,aa)+square_covariance(v,v,bb)
              -square_covariance(u,v,ab)-square_covariance(v,u,ab.T))
    return expected,(variance+variance.T)/2


def equal_weights(sp):
    ids=np.unique(sp);w=np.zeros(len(sp))
    for s in ids:
        mask=sp==s;w[mask]=1/(len(ids)*mask.sum())
    return w


def component_sizes(sp,st):
    groups={s:set(sp[st==s]) for s in np.unique(st)}
    pending=set(sp);sizes=[]
    while pending:
        reached={min(pending)}
        while True:
            nxt=set(reached)
            for members in groups.values():
                if reached & members: nxt.update(members)
            if nxt==reached: break
            reached=nxt
        pending-=reached;sizes.append(len(reached))
    return sorted(sizes,reverse=True)


def analyse(X,H,sp,st,role):
    sp=np.asarray(sp);st=np.asarray(st);train=np.asarray(role)=='train'
    if X.shape!=(165,21) or H.shape!=(165,) or (train.sum(),(~train).sum())!=(91,74):
        raise ValueError('frozen geometry dimensions changed')
    if (len(np.unique(sp[train])),len(np.unique(sp[~train])))!=(29,27):
        raise ValueError('frozen species counts changed')
    if set(st[train])&set(st[~train]) or set(sp[train])&set(sp[~train]):
        raise ValueError('original study/species holdout leakage')
    w=equal_weights(sp[train]);wv=equal_weights(sp[~train])
    P=ridge_operator(X,train,w);Q=ridge_operator(np.column_stack([X,H]),train,w)
    A=error_operator(P,train);B=error_operator(Q,train)
    mu,sd=scale(X[train],w);Z=(X-mu)/sd
    coeff=np.array([(0.12 if i%2==0 else -0.12) if i<11 else
                    (0.06 if i%2==0 else -0.06) for i in range(21)])
    m=Z@coeff;omega=noise_covariance(sp,st)
    bm=A@m;hm=B@m
    bvar=np.diag(A@omega@A.T);hvar=np.diag(B@omega@B.T)
    loss_b=float(wv@(bm*bm+bvar));loss_h=float(wv@(hm*hm+hvar))
    expected,covdiff=loss_difference_moments(A,B,m,omega)
    ids=np.unique(sp[~train]);G=np.zeros((len(ids),(~train).sum()))
    for i,s in enumerate(ids):
        mask=sp[~train]==s;G[i,mask]=1/mask.sum()
    C=G@covdiff@G.T
    fullvar=float(C.sum()/len(ids)**2);diagvar=float(np.trace(C)/len(ids)**2)
    heldsp=sp[~train];heldst=st[~train]
    shared=[s for s in np.unique(heldst) if len(set(heldsp[heldst==s]))>1]
    shared_pairs=set()
    for study in shared:
        members=sorted(set(heldsp[heldst==study]))
        for i,s in enumerate(members):
            for t in members[i+1:]: shared_pairs.add((s,t))
    meantrain_cov=np.zeros_like(omega)
    meantrain_cov[np.ix_(~train,~train)]=omega[np.ix_(~train,~train)]
    _,conditional=loss_difference_moments(A,B,m,meantrain_cov)
    C0=G@conditional@G.T
    conditional_full=float(C0.sum()/len(ids)**2)
    conditional_diag=float(np.trace(C0)/len(ids)**2)
    return {
        'source_train_species':29,'source_heldout_species':27,
        'heldout_source_studies':len(np.unique(heldst)),
        'heldout_studies_containing_multiple_species':len(shared),
        'heldout_species_pairs_sharing_a_study':len(shared_pairs),
        'heldout_study_species_component_sizes':component_sizes(heldsp,heldst),
        'train_test_study_overlap':0,
        'analytic_null_expected_baseline_MSE':loss_b,
        'analytic_null_expected_augmented_MSE':loss_h,
        'ratio_of_expected_MSE_reduction_NOT_expected_ratio':(loss_b-loss_h)/loss_b,
        'baseline_bias_squared':float(wv@(bm*bm)),
        'augmented_bias_squared':float(wv@(hm*hm)),
        'baseline_prediction_error_variance':float(wv@bvar),
        'augmented_prediction_error_variance':float(wv@hvar),
        'mean_loss_difference_integrating_training_and_test':float(wv@expected),
        'unconditional_loss_difference_variance':fullvar,
        'unconditional_diagonal_species_variance':diagvar,
        'unconditional_full_to_diagonal_variance_ratio':fullvar/diagvar if diagvar>0 else None,
        'heldout_only_covariance_at_mean_training_full_variance':conditional_full,
        'heldout_only_covariance_at_mean_training_diagonal_variance':conditional_diag,
        'heldout_only_full_to_diagonal_variance_ratio':conditional_full/conditional_diag if conditional_diag>0 else None,
        'variance_ratios_are_NOT_bootstrap_coverage_estimates':True,
        'cannot_attribute_seven_saved_decisions_without_trial_level_outputs':True,
    }


def source_reconstruction(work):
    # Import source-only routines, NEVER known_truth_gate/run or empirical scripts.
    from proceed70_haldane_source_null_calibration import source_only_geometry
    from proceed_v62_borealis_header_only import validated_get_catalog,verify_frozen_catalog
    from proceed_v62_metadata_row_census import fetch_original_datafile_bytes
    from proceed_v62_genetic_geometry_intersection import SOURCE_CONTRACT,NATIVE_LOCK
    from proceed_plant_chelsa_source_support import source_geometry,verified_raster
    from site_feature_core import raster_values
    here=Path(__file__).resolve().parent
    native=json.loads(NATIVE_LOCK.read_text())['native_original_transport']
    catalog=validated_get_catalog('https://borealisdata.ca/api/datasets/:persistentId/?persistentId=doi%3A10.5683%2FSP3%2FNXSL3Q')
    verify_frozen_catalog(catalog,json.loads(SOURCE_CONTRACT.read_text())['current_source'])
    raw,md5,sha=fetch_original_datafile_bytes(1024206)
    if len(raw)!=native['observed_size_bytes'] or md5!=native['observed_md5'] or sha!=native['observed_sha256']:
        raise ValueError('original source checksum failed')
    design=json.loads((here/'proceed_plant_chelsa_source_support_contract_v01.json').read_text())
    full,_=source_geometry(raw,design)
    coords=sorted({v for _,_,a,b in full for v in (a,b)})
    xy=np.asarray(coords,float)
    digest=lambda a:hashlib.sha256(np.asarray(a,dtype='<f8').tobytes()).hexdigest()
    if digest(xy)!=FROZEN_COORD:raise ValueError('source coordinates changed')
    assets=json.loads((here/'chelsa_trace21k_fixed8_assets_v01.json').read_text())['assets']
    columns=[]
    for i,asset in enumerate(assets):
        p=work/f'audit_{i}.tif';verified_raster(asset,p)
        columns.append(raster_values(p,xy));p.unlink()
        print('UNCHANGED_SOURCE',i+1,flush=True)
    matrix=np.column_stack(columns)
    if digest(matrix)!=FROZEN_CLIMATE:raise ValueError('source climate matrix changed')
    X,H,sp,st,role,_,_=source_only_geometry(raw,matrix,coords,design)
    if digest(X)!=FROZEN_X or digest(H)!=FROZEN_H:raise ValueError('frozen source model matrices changed')
    return X,H,sp,st,role


def main():
    p=argparse.ArgumentParser();p.add_argument('--work',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True);args=p.parse_args()
    here=Path(__file__).resolve().parent
    terminal=json.loads((here/'proceed70_haldane_synthetic_terminal_receipt_v01.json').read_text())
    status=terminal['status']
    if status!='TERMINAL_NOT_EVALUABLE_H_SYNTHETIC_NULL_CEILING_FAIL_NO_BIOLOGICAL_OUTCOME_OPENING':
        raise ValueError('this diagnostic is authorized ONLY after terminal closure')
    report={'schema':'ttf_terminal_haldane_closed_form_diagnosis_v01',
            'terminal_status_unchanged':status,'real_phenotype_rate_values_read':False,
            'additional_null_or_positive_worlds_drawn':0,
            'parameters_or_decision_thresholds_changed':False,
            'diagnostic_not_new_calibration_or_authorization':True}
    args.work.mkdir(parents=True,exist_ok=True);args.output.parent.mkdir(parents=True,exist_ok=True)
    try:
        report.update(analyse(*source_reconstruction(args.work)))
        report['saved_7_of_100_Wilson95_recomputed']=wilson(7,100)
        report['max_allowed_positives_among_100_for_original_Wilson_rule']=max(k for k in range(101) if wilson(k,100)[1]<=.1)
        report['duplicate_covariate_algebra_gain_without_new_information']=1-((1+ALPHA)/(2+ALPHA))**2
        report['status']='ANALYTICAL_DIAGNOSIS_COMPLETE_TERMINAL_STOP_UNCHANGED'
    except Exception as exc:
        report.update(status='DIAGNOSTIC_ERROR_TERMINAL_STOP_UNCHANGED',error=f'{type(exc).__name__}: {exc}')
        args.output.write_text(json.dumps(report,indent=2,sort_keys=True)+'\n')
        raise
    args.output.write_text(json.dumps(report,indent=2,sort_keys=True)+'\n')
    print(json.dumps(report,sort_keys=True))

if __name__=='__main__':main()
