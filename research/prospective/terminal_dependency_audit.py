"""Post-terminal source dependency inventory. No climate rebuild or new worlds.

The original Gaussian null's covariance is integrated only for species MEAN
labels, not squared prediction errors. This does not estimate bootstrap
coverage or explain any particular saved rejection.
"""
from __future__ import annotations
import argparse,json
from collections import defaultdict,Counter
from pathlib import Path


def summarize(species,studies,roles):
    if not(len(species)==len(studies)==len(roles)):
        raise ValueError('unaligned source metadata')
    trs={s for s,r in zip(species,roles) if r=='train'}
    tes={s for s,r in zip(species,roles) if r=='heldout'}
    trst={t for t,r in zip(studies,roles) if r=='train'}
    test={t for t,r in zip(studies,roles) if r=='heldout'}
    if trs&tes or trst&test:raise ValueError('study/species holdout leakage')
    by_study=defaultdict(set); by_species=defaultdict(Counter)
    for s,t,r in zip(species,studies,roles):
        if r=='heldout':by_study[t].add(s);by_species[s][t]+=1
    shared=[v for v in by_study.values() if len(v)>1]
    pairs=set()
    for members in shared:
        ss=sorted(members)
        pairs.update((a,b) for i,a in enumerate(ss) for b in ss[i+1:])
    pending=set(tes);components=[]
    while pending:
        members={min(pending)}
        while True:
            extra=set(members)
            for g in by_study.values():
                if g&members:extra.update(g)
            if extra==members:break
            members=extra
        pending-=members;components.append(len(members))
    # Mean of original synthetic labels for each species:
    # one .25 species intercept, .25 study effects and .25 independent noise.
    diag=0.;off=0.
    for s,count in by_species.items():
        n=sum(count.values())
        diag+=.25+.25*sum((v/n)**2 for v in count.values())+.25/n
    for s,t in pairs:
        c=by_species[s];d=by_species[t];ns=sum(c.values());nt=sum(d.values())
        off+=2*.25*sum(c[k]*d.get(k,0)/(ns*nt) for k in c)
    return dict(train_species=len(trs),heldout_species=len(tes),
        train_studies=len(trst),heldout_studies=len(test),
        train_test_species_overlap=0,train_test_study_overlap=0,
        heldout_studies_with_multiple_species=len(shared),
        heldout_species_pairs_sharing_a_study=len(pairs),
        heldout_connected_component_sizes=sorted(components,reverse=True),
        synthetic_heldout_species_mean_covariance_diagonal_sum=diag,
        synthetic_heldout_species_mean_covariance_offdiagonal_sum=off,
        synthetic_species_mean_average_variance=(diag+off)/len(tes)**2 if tes else None,
        synthetic_diagonal_only_mean_average_variance=diag/len(tes)**2 if tes else None,
        covariance_is_for_synthetic_species_mean_labels_NOT_prediction_loss=True,
        bootstrap_error_rate_or_causal_attribution_estimated=False)


def main():
    from proceed70_haldane_source_null_calibration import presence_keyed_metadata,PRESENCE
    from proceed_v62_metadata_row_census import fetch_original_datafile_bytes
    from proceed_plant_chelsa_source_support import source_geometry,component_assignment
    from proceed_v62_genetic_geometry_intersection import NATIVE_LOCK
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True)
    args=p.parse_args();here=Path(__file__).resolve().parent
    terminal=json.loads((here/'proceed70_haldane_synthetic_terminal_receipt_v01.json').read_text())
    if terminal['status']!='TERMINAL_NOT_EVALUABLE_H_SYNTHETIC_NULL_CEILING_FAIL_NO_BIOLOGICAL_OUTCOME_OPENING':
        raise ValueError('terminal status must remain closed')
    lock=json.loads(NATIVE_LOCK.read_text())['native_original_transport']
    raw,md5,sha=fetch_original_datafile_bytes(1024206)
    if md5!=lock['observed_md5'] or sha!=lock['observed_sha256'] or len(raw)!=lock['observed_size_bytes']:
        raise ValueError('original source identity mismatch')
    full,graph=source_geometry(raw,{})
    mapping,split,digest=component_assignment(graph)
    if digest!='73b16a46f4458aee64d735a495273e5d760683b020c5432a9701fc14214f231d':
        raise ValueError('original source roles changed')
    meta=presence_keyed_metadata(raw,json.loads(PRESENCE.read_text()))
    keys=sorted(meta)
    if len(full)!=1561 or len(keys)!=165:raise ValueError('source row geometry changed')
    sp=[k[0] for k in keys];st=[k[1] for k in keys];role=[mapping[s] for s in sp]
    out=summarize(sp,st,role)
    if (out['train_species'],out['heldout_species'])!=(29,27):raise ValueError('masked taxa changed')
    out.update(schema='ttf_terminal_haldane_source_dependency_only_v01',
        status='METADATA_DEPENDENCY_DIAGNOSIS_COMPLETE_TERMINAL_STOP_UNCHANGED',
        original_source_sha256=sha,original_source_role_sha256=digest,
        original_mask_pair_count=len(keys),additional_worlds_drawn=0,
        numeric_phenotype_or_Haldane_values_parsed=False,
        changed_alpha_or_null_threshold=False,
        original_numeric_model_risk_not_computed=True,
        original_model_matrix_identity_issue='Run 43594694027 recovered original 8-band climate bytes but failed 21-control/history matrix hashes. Exact original-model expected-risk calculation remains unavailable; no acceptance of new hashes.',
        original_terminal_decision_unchanged=terminal['status'])
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(out,sort_keys=True,indent=2)+'\n');print(json.dumps(out,sort_keys=True))

if __name__=='__main__':main()
