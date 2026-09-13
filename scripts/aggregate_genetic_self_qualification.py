#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
from pathlib import Path

from ttf.precision import wilson_interval


def main() -> int:
    ap=argparse.ArgumentParser()
    ap.add_argument('--input-dir',type=Path,required=True)
    ap.add_argument('--qualification',type=Path,required=True)
    ap.add_argument('--output',type=Path,required=True)
    args=ap.parse_args()
    qual=json.loads(args.qualification.read_text())
    expected={'null':int(qual['synthetic_worlds']['null_evaluation_worlds']),
              'private_A2':int(qual['synthetic_worlds']['private_A2_evaluation_worlds'])}
    grouped={cell:{} for cell in expected}
    for path in sorted(args.input_dir.rglob('*.json')):
        p=json.loads(path.read_text())
        if p.get('schema')!='ttf_genetic_self_evaluation_shard_v0.1':
            continue
        if p.get('empirical_genetic_outcomes_opened') is not False:
            raise RuntimeError('empirical outcome firewall drift')
        cell=str(p['cell'])
        if cell not in grouped:
            raise RuntimeError(f'unexpected self cell {cell}')
        for row in p['rows']:
            rep=int(row['replicate'])
            if rep in grouped[cell]:
                raise RuntimeError(f'duplicate self replicate {cell}:{rep}')
            grouped[cell][rep]=float(row['p_value'])
    alpha=float(qual['inference']['alpha'])
    cells={}
    for cell,n in expected.items():
        if sorted(grouped[cell])!=list(range(n)):
            missing=sorted(set(range(n))-set(grouped[cell]))[:10]
            raise RuntimeError(f'incomplete self cell {cell}; first missing={missing}')
        reject=sum(grouped[cell][i] <= alpha for i in range(n))
        interval=wilson_interval(reject,n)
        cells[cell]={
            'n_worlds':n,'rejections':int(reject),'rejection_rate':float(reject/n),
            'wilson95_low':float(interval.low),'wilson95_high':float(interval.high),
        }
    type1_pass=cells['null']['wilson95_high'] <= 0.10
    power_pass=cells['private_A2']['wilson95_low'] >= 0.80
    passed=bool(type1_pass and power_pass)
    out={
        'schema':'ttf_genetic_self_detectability_qualification_result_v0.1',
        'status':'PASS' if passed else 'SELF_DETECTABILITY_NOT_QUALIFIED',
        'cells':cells,
        'type1_gate':{'wilson95_upper_ceiling':0.10,'pass':bool(type1_pass)},
        'power_gate':{'wilson95_lower_floor':0.80,'pass':bool(power_pass)},
        'passed':passed,
        'empirical_genetic_outcomes_opened':False,
        'claim_boundary':'PASS only qualifies the within-species detectability diagnostic. It does not open or analyze empirical genetic outcomes.',
    }
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(out,indent=2,sort_keys=True)+'\n')
    print(json.dumps(out,sort_keys=True))
    return 0

if __name__=='__main__':
    raise SystemExit(main())
