#!/usr/bin/env python3
from __future__ import annotations

import argparse,json
from pathlib import Path

from ttf.precision import wilson_interval


def collect(directory:Path,panel:str,cell:str,expected:int):
    rows=[]
    for path in sorted(directory.glob("*.json")):
        p=json.loads(path.read_text())
        if p.get("schema")!="ttf_lepidoptera_postprimary_self_shard_v0.1":
            continue
        if p.get("panel")!=panel or p.get("cell")!=cell:
            continue
        start=int(p["start"])
        stats=list(map(float,p["statistics"]))
        if len(stats)!=int(p["count"]):
            raise RuntimeError("self shard count drift")
        pvals=p.get("p_values")
        if cell=="reference":
            rows.extend((start+i,stats[i],None) for i in range(len(stats)))
        else:
            if pvals is None or len(pvals)!=len(stats):
                raise RuntimeError("self evaluation p-value drift")
            rows.extend((start+i,stats[i],float(pvals[i])) for i in range(len(stats)))
    if len(rows)!=expected:
        raise RuntimeError(f"{panel}/{cell}: {len(rows)} != {expected}")
    rows.sort()
    if [i for i,_,_ in rows]!=list(range(expected)):
        raise RuntimeError(f"{panel}/{cell}: replicate coverage drift")
    return rows


def main()->int:
    ap=argparse.ArgumentParser()
    ap.add_argument("--input-dir",type=Path,required=True)
    ap.add_argument("--rule",type=Path,required=True)
    ap.add_argument("--panel",choices=("butterfly","host"),required=True)
    ap.add_argument("--output",type=Path,required=True)
    args=ap.parse_args()

    rule=json.loads(args.rule.read_text())
    if rule.get("schema")!="ttf_lepidoptera_postprimary_self_detectability_rule_v0.1":
        raise RuntimeError("self rule schema drift")
    ref=collect(args.input_dir,args.panel,"reference",1000)
    null=collect(args.input_dir,args.panel,"null",500)
    pos=collect(args.input_dir,args.panel,"private_A2",500)

    alpha=float(rule["method_inheritance"]["alpha"])
    def cell(rows):
        p=[r[2] for r in rows]
        k=sum(v<=alpha for v in p)
        ci=wilson_interval(k,len(rows))
        stats=[r[1] for r in rows]
        return {
            "worlds":len(rows),
            "rejections":int(k),
            "rejection_rate":float(k/len(rows)),
            "wilson95_lower":float(ci.low),
            "wilson95_upper":float(ci.high),
            "mean_statistic":float(sum(stats)/len(stats)),
        }

    null_cell=cell(null)
    pos_cell=cell(pos)
    gate=rule["exact_geometry_qualification_per_panel"]
    type1_pass=null_cell["wilson95_upper"]<=0.10
    power_pass=pos_cell["wilson95_lower"]>=0.80
    passed=bool(type1_pass and power_pass)
    ref_stats=[r[1] for r in ref]

    out={
        "schema":"ttf_lepidoptera_postprimary_self_qualification_v0.1",
        "status":"PASS" if passed else "SELF_DETECTABILITY_NOT_QUALIFIED",
        "panel":args.panel,
        "reference":{
            "worlds":1000,
            "mean_statistic":float(sum(ref_stats)/len(ref_stats)),
        },
        "evaluation":{
            "null":null_cell,
            "private_A2":pos_cell,
        },
        "gates":{
            "type1_pass":bool(type1_pass),
            "private_A2_power_pass":bool(power_pass),
            "all_pass":passed,
        },
        "empirical_self_eligible":passed,
        "cross_species_primary_rewritten":False,
    }
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(out,indent=2,sort_keys=True)+"\n")
    print(json.dumps(out["gates"],sort_keys=True))
    return 0


if __name__=="__main__":
    raise SystemExit(main())
