#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
from pathlib import Path


def main() -> int:
    ap=argparse.ArgumentParser()
    ap.add_argument('--input-dir',type=Path,required=True)
    ap.add_argument('--qualification',type=Path,required=True)
    ap.add_argument('--output',type=Path,required=True)
    args=ap.parse_args()
    qual=json.loads(args.qualification.read_text())
    expected=int(qual['synthetic_worlds']['independent_null_reference_worlds'])
    values={}
    for path in sorted(args.input_dir.rglob('*.json')):
        p=json.loads(path.read_text())
        if p.get('schema')!='ttf_genetic_self_reference_shard_v0.1':
            continue
        if p.get('empirical_genetic_outcomes_opened') is not False:
            raise RuntimeError('empirical outcome firewall drift')
        for row in p['rows']:
            rep=int(row['replicate'])
            if rep in values:
                raise RuntimeError(f'duplicate self reference replicate {rep}')
            values[rep]=float(row['statistic'])
    if sorted(values)!=list(range(expected)):
        missing=sorted(set(range(expected))-set(values))[:10]
        raise RuntimeError(f'incomplete self reference; first missing={missing}')
    out={
        'schema':'ttf_genetic_self_references_v0.1',
        'status':'independent_null_reference_complete',
        'n_worlds':expected,
        'statistics':[values[i] for i in range(expected)],
        'empirical_genetic_outcomes_opened':False,
        'qualification_claim_made':False,
    }
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(out,indent=2,sort_keys=True)+'\n')
    print(json.dumps({k:v for k,v in out.items() if k!='statistics'},sort_keys=True))
    return 0

if __name__=='__main__':
    raise SystemExit(main())
