#!/usr/bin/env python3
from __future__ import annotations

import argparse,json
from collections import Counter
from pathlib import Path
import numpy as np

from ttf.calibration import CalibrationCell
from ttf.precision import qualify_calibration_precision
from ttf.profiled_private_null import profiled_private_pvalue

LABELS=('A0','A0p5','A1','A2','A3','A5','A10','infinite_snr')
MANDATORY=((0.0,0.5),(0.0,1.0),(0.0,2.0),(0.0,3.0),(1.0,2.0))
MATCH={0.5:'A0p5',1.0:'A1',2.0:'A2',3.0:'A3'}


def main()->int:
    ap=argparse.ArgumentParser(); ap.add_argument('--input-dir',type=Path,required=True); ap.add_argument('--rule',type=Path,required=True); ap.add_argument('--source-ledger',type=Path,required=True); ap.add_argument('--output',type=Path,required=True); args=ap.parse_args()
    rule=json.loads(args.rule.read_text()); source=json.loads(args.source_ledger.read_text())
    if rule.get('status')!='frozen_before_v08_external_validity_panels_are_declared': raise RuntimeError('v08 architecture drift')
    if source.get('synthetic_worlds_run_at_freeze')!=0 or source.get('empirical_trait_or_colour_fields_read') is not False: raise RuntimeError('reserve not unopened')
    refs={}; obs={}
    for path in sorted(args.input_dir.glob('*.json')):
        d=json.loads(path.read_text())
        if d.get('schema')!='ttf_v08_gate_d_statistics_v0.1': continue
        if d.get('empirical_colour_outcome_opened') is not False: raise RuntimeError('empirical colour firewall failed')
        cfg=d['config']
        if d['mode']=='reference':
            label=cfg['configuration_label'];
            if label in refs: raise RuntimeError('duplicate ref')
            if len(d['statistics'])!=1999: raise RuntimeError('ref size drift')
            refs[label]=d
        else:
            key=(float(cfg['shared_fraction']),float(cfg['amplitude']))
            if key in obs: raise RuntimeError('duplicate observed')
            if len(d['statistics'])!=500: raise RuntimeError('observed size drift')
            obs[key]=d
    if set(refs)!=set(LABELS): raise RuntimeError(('reference set',sorted(refs)))
    if set(obs)!=set(MANDATORY): raise RuntimeError(('observed set',sorted(obs)))
    reference={label:(np.asarray(d['training_strength'],float),np.asarray(d['statistics'],float)) for label,d in refs.items()}
    cells=[]; recovery={}; pair_counts={}
    for shared,amp in MANDATORY:
        d=obs[(shared,amp)]; t=np.asarray(d['statistics'],float); u=np.asarray(d['training_strength'],float); p=[]; pairs=[]
        for ti,ui in zip(t,u):
            pv,selected,_,_=profiled_private_pvalue(float(ti),float(ui),reference,profile_draws=999,selected_configs=2,scale_floor=1e-6)
            p.append(pv); pairs.append('+'.join(selected))
        p=np.asarray(p,float)
        cells.append(CalibrationCell(shared_fraction=shared,amplitude=amp,n_replicates=500,alpha=.05,rejection_rate=float(np.mean(p<=.05)),mean_statistic=float(t.mean()),mean_null_statistic=float('nan'),median_p_value=float(np.median(p))))
        key=f's{shared}-a{amp}'; pair_counts[key]={k:int(v) for k,v in sorted(Counter(pairs).items())}
        recovery[key]=None if shared>0 else float(np.mean([MATCH[amp] in x.split('+') for x in pairs]))
    precision=qualify_calibration_precision(cells,moderate_amplitude=2.0,type1_upper_ceiling=.10,power_lower_floor=.80).to_dict()
    out={'schema':'ttf_v08_gate_d_rgfca_reserve_v0.1','status':'synthetic_deployment_adequacy_on_unopened_rgfca_reserve','rule':str(args.rule),'source_ledger':str(args.source_ledger),'cells':[c.to_dict() for c in cells],'precision_qualification':precision,'gate_d_pass':bool(precision['passed']),'matching_private_label_recovery':recovery,'selected_pair_counts':pair_counts,'empirical_colour_outcome_opened':False,'claim_ready':False,'interpretation_boundary':'Gate-D is a pre-outcome design-adequacy check. A pass permits subsequent empirical TTF analysis on this frozen reserve design; it does not itself constitute an empirical flower-colour result.'}
    args.output.parent.mkdir(parents=True,exist_ok=True); args.output.write_text(json.dumps(out,indent=2,sort_keys=True)+'\n'); print(json.dumps({'gate_d_pass':out['gate_d_pass'],'precision':precision},sort_keys=True)); return 0

if __name__=='__main__': raise SystemExit(main())
