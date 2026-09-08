#!/usr/bin/env python3
from __future__ import annotations

import argparse,json
from collections import Counter
from pathlib import Path
import numpy as np

from ttf.calibration import CalibrationCell
from ttf.precision import wilson_interval
from ttf.profiled_private_null import profiled_private_pvalue

PANELS=('europe_birds','na_waterbirds')
PRIVATE_AMPLITUDES=(0.5,1.0,2.0,3.0)
LABELS=('A0','A0p5','A1','A2','A3','A5','A10','infinite_snr')
MATCH={0.5:'A0p5',1.0:'A1',2.0:'A2',3.0:'A3'}


def main()->int:
    ap=argparse.ArgumentParser()
    ap.add_argument('--input-dir',type=Path,required=True)
    ap.add_argument('--rule',type=Path,required=True)
    ap.add_argument('--panel-source',type=Path,required=True)
    ap.add_argument('--output',type=Path,required=True)
    args=ap.parse_args()
    rule=json.loads(args.rule.read_text())
    source=json.loads(args.panel_source.read_text())
    if rule.get('status')!='frozen_before_v08_external_validity_panels_are_declared': raise RuntimeError('v08 rule drift')
    if source.get('status')!='two_fresh_external_validity_panels_frozen_after_v08_rule': raise RuntimeError('v08 panels not prospectively frozen')
    if source.get('synthetic_worlds_run_at_freeze')!=0 or source.get('candidate_performance_evaluated_at_freeze') is not False: raise RuntimeError('panel freeze firewall failed')
    if rule['gate_v_external_validity']['shared_positive_control_forbidden_in_gate_v'] is not True: raise RuntimeError('shared positive firewall missing')

    payloads=[]
    for path in sorted(args.input_dir.glob('*.json')):
        d=json.loads(path.read_text())
        if d.get('schema')=='ttf_v08_validity_statistics_v0.1': payloads.append(d)
    reference={p:{} for p in PANELS}; observed={p:{} for p in PANELS}
    for d in payloads:
        if d.get('shared_positive_control_opened') is not False or d.get('legacy_birds_butterflies_opened') is not False or d.get('rgfca_reserve_opened') is not False: raise RuntimeError('v08 validity firewall failure')
        p=d['panel']; cfg=d['config']
        if p not in PANELS or float(cfg['shared_fraction'])!=0.0: raise RuntimeError('unexpected panel/shared state')
        if d['mode']=='reference':
            label=cfg['configuration_label']
            if label in reference[p]: raise RuntimeError('duplicate reference')
            if len(d['statistics'])!=1999 or len(d['training_strength'])!=1999: raise RuntimeError('reference size drift')
            reference[p][label]=d
        else:
            a=float(cfg['amplitude'])
            if a in observed[p]: raise RuntimeError('duplicate observed')
            if len(d['statistics'])!=500 or len(d['training_strength'])!=500: raise RuntimeError('observed size drift')
            observed[p][a]=d
    for p in PANELS:
        if set(reference[p])!=set(LABELS): raise RuntimeError((p,'reference set',sorted(reference[p])))
        if set(observed[p])!=set(PRIVATE_AMPLITUDES): raise RuntimeError((p,'observed set',sorted(observed[p])))

    reports={}; both=True
    for p in PANELS:
        refs={label:(np.asarray(d['training_strength'],float),np.asarray(d['statistics'],float)) for label,d in reference[p].items()}
        cells=[]; recovery={}; selected_counts={}
        for a in PRIVATE_AMPLITUDES:
            d=observed[p][a]; t=np.asarray(d['statistics'],float); u=np.asarray(d['training_strength'],float)
            ps=[]; pairs=[]
            for ti,ui in zip(t,u):
                pv,selected,_,_=profiled_private_pvalue(float(ti),float(ui),refs,profile_draws=999,selected_configs=2,scale_floor=1e-6)
                ps.append(pv); pairs.append('+'.join(selected))
            ps=np.asarray(ps,float)
            cells.append(CalibrationCell(shared_fraction=0.0,amplitude=a,n_replicates=500,alpha=.05,rejection_rate=float(np.mean(ps<=.05)),mean_statistic=float(t.mean()),mean_null_statistic=float('nan'),median_p_value=float(np.median(ps))))
            matching=MATCH[a]
            recovery[str(a)]=float(np.mean([matching in pair.split('+') for pair in pairs]))
            selected_counts[str(a)]={k:int(v) for k,v in sorted(Counter(pairs).items())}
        intervals=[]
        for c in cells:
            successes=int(round(c.rejection_rate*c.n_replicates))
            interval=wilson_interval(successes,c.n_replicates)
            intervals.append({'amplitude':c.amplitude,'estimate':c.rejection_rate,'low':interval.low,'high':interval.high,'successes':successes,'trials':c.n_replicates})
        max_high=max(x['high'] for x in intervals); max_est=max(x['estimate'] for x in intervals)
        passed=bool(max_high<=0.10); both=both and passed
        reports[p]={'cells':[c.to_dict() for c in cells],'zero_shared_intervals':intervals,'max_zero_shared_estimate':max_est,'max_zero_shared_upper95':max_high,'type1_upper_ceiling':0.10,'passed':passed,'matching_private_label_recovery':recovery,'selected_pair_counts':selected_counts,'shared_positive_control_opened':False}

    out={'schema':'ttf_v08_external_validity_v0.1','status':'fresh_external_type1_validity_confirmation','rule':str(args.rule),'panel_source':str(args.panel_source),'panels':reports,'gate_v_pass':bool(both),'selection_rule_applied_mechanically':True,'power_claim_made':False,'shared_positive_control_opened':False,'legacy_birds_butterflies_opened':False,'rgfca_reserve_opened':False,'next_gate':'If gate_v_pass is true, the RGFCA reserve may receive the separately frozen Gate-D local type-I plus shared-A2 detectability audit before any empirical colour outcome is opened. If false, preserve failure and do not open the reserve.','claim_ready':False}
    args.output.parent.mkdir(parents=True,exist_ok=True); args.output.write_text(json.dumps(out,indent=2,sort_keys=True)+'\n')
    print(json.dumps({'gate_v_pass':both,'panels':{p:{'max_upper95':reports[p]['max_zero_shared_upper95'],'max_estimate':reports[p]['max_zero_shared_estimate'],'recovery':reports[p]['matching_private_label_recovery']} for p in PANELS}},sort_keys=True))
    return 0

if __name__=='__main__': raise SystemExit(main())
