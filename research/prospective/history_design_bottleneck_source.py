"""Anatomy of the frozen source design; no empirical rate or new null worlds.

Continue only on BOTH original X/H hashes. Report every declared nuisance
span and idealized new-block planning; no method is chosen to reopen v1.
"""
from __future__ import annotations
import argparse
import hashlib
import json
from pathlib import Path
import tempfile
import numpy as np
from history_design_bottleneck import analyse_geometry, required_blocks, ALPHA

HERE=Path(__file__).resolve().parent
CLOSED='TERMINAL_NOT_EVALUABLE_H_SYNTHETIC_NULL_CEILING_FAIL_NO_BIOLOGICAL_OUTCOME_OPENING'


def digest(a):
    return hashlib.sha256(np.asarray(a,dtype='<f8').tobytes(order='C')).hexdigest()


def require_identity(X,H,xsha,hsha):
    if (X.shape!=(165,21) or H.shape!=(165,) or
            not np.isfinite(X).all() or not np.isfinite(H).all() or
            digest(X)!=xsha or digest(H)!=hsha):
        raise ValueError('Original input shape, finiteness or BOTH exact byte hashes differ')
    return True


def inspect_source(X,H,sp,st,role):
    from history_absorbed_source import analyse_source
    # This existing check also enforces original study/species split and shape.
    prior=analyse_source(X,H,sp,st,role)
    train=role=='train'; held=role=='heldout'
    w=np.zeros(165)
    for s in np.unique(sp[train]):
        take=train&(sp==s);w[take]=1/(29*int(take.sum()))
    hm=float(w@H);hs=float(np.sqrt(w@((H-hm)**2)))
    xm=w@X;xs=np.sqrt(w@((X-xm)**2));xs=np.where(xs>1e-10,xs,1.)
    r=analyse_geometry((X[held]-xm)/xs,(H[held]-hm)/hs,sp[held],st[held])
    final=r['spans']['groups_baseline21'];old=prior['within_group_design']
    if (final['rank']!=old['nuisance_rank'] or
            final['residual_df_after_history']!=old['residual_df'] or
            not np.isclose(final['history_energy'],old['history_residual_energy'],rtol=1e-8,atol=1e-11)):
        raise ValueError('New diagnostic does not reproduce prior final span')
    r['previous_final_design_reproduced']=True
    r['original_training_history_SD']=hs
    return r


def main():
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True)
    args=p.parse_args();args.output.parent.mkdir(parents=True,exist_ok=True)
    terminal_path=HERE/'proceed70_haldane_synthetic_terminal_receipt_v01.json'
    before=terminal_path.read_bytes();terminal=json.loads(before)
    if terminal['status']!=CLOSED:raise ValueError('Old biological STOP must remain')
    out={'schema':'ttf_history_design_bottleneck_and_new_contrast_targets_v01',
         'actual_trait_or_fitness_outcomes_read':False,'new_random_worlds':0,
         'empirical_opening_authorized':False,'old_terminal_status':CLOSED,
         'method_selected_from_comparison':False}
    try:
        from terminal_numeric_provenance import code_identity,prepare_sources
        from history_absorbed_source import bounded_runtime_reconstruction,select_exact_mode
        out['code_identity']=code_identity()
        if not all(v['identical'] for v in out['code_identity'].values()):
            raise ValueError('Original source helpers modified')
        with tempfile.TemporaryDirectory(prefix='ttf-private-anatomy-') as td:
            work=Path(td);out['source_identity']=prepare_sources(work,terminal)
            runtime={};bounded_runtime_reconstruction(work,runtime)
            lock=terminal['source_pre_registered']
            xsha=lock['source_21_baseline_numeric_SHA256'];hsha=lock['source_one_history_numeric_SHA256']
            selected=select_exact_mode(runtime['runtime_modes'],xsha,hsha)
            out['selected_exact_runtime']=selected
            out['runtime_hashes']={k:{q:v.get(q) for q in ('X_sha256','H_sha256','worker_failed')}
                                   for k,v in runtime['runtime_modes'].items()}
            if selected is None:raise ValueError('No bounded runtime reproduces both original hashes')
            with np.load(work/(selected+'.npz'),allow_pickle=False) as a:
                require_identity(a['X'],a['H'],xsha,hsha)
                out['design_anatomy']=inspect_source(a['X'],a['H'],a['sp'],a['st'],a['role'])
            out['source_X_sha256']=xsha;out['source_H_sha256']=hsha
        out['idealized_new_independent_design']={
            'alpha':ALPHA,'power_target':.8,'gamma':.5,'iid_normal_sigma':.5,
            'two_observations_per_block':True,
            'baseline_covariates_identical_within_each_new_block':True,
            'independent_errors_and_common_history_coefficient_assumed':True,
            'current_data_or_additional_Haldane_results_not_used':True,
            'rows':[required_blocks(d) for d in (.5,1.,2.)],
            'not_a_real_world_guaranteed_sample_size':True}
        if before!=terminal_path.read_bytes():raise ValueError('Old terminal decision was changed')
        out['status']='DESIGN_BOTTLENECK_DIAGNOSED_OLD_RESPONSE_STOP_UNCHANGED'
    except Exception as exc:
        out['status']='SOURCE_OR_DESIGN_DIAGNOSTIC_STOP_OLD_RESPONSE_STOP_UNCHANGED'
        out['error']=type(exc).__name__+': '+str(exc)[:450]
        args.output.write_text(json.dumps(out,indent=2,sort_keys=True)+'\n')
        raise
    args.output.write_text(json.dumps(out,indent=2,sort_keys=True)+'\n')
    print(json.dumps(out,sort_keys=True))

if __name__=='__main__':main()
