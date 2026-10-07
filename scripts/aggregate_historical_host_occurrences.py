#!/usr/bin/env python3
from __future__ import annotations

import argparse,csv,hashlib,json,math
from pathlib import Path


def sha256_path(path: Path) -> str:
    h=hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda:f.read(1<<20),b""): h.update(chunk)
    return h.hexdigest()


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--panel",type=Path,required=True)
    ap.add_argument("--rule",type=Path,required=True)
    ap.add_argument("--shard-dir",type=Path,required=True)
    ap.add_argument("--output-dir",type=Path,required=True)
    args=ap.parse_args()
    rule=json.loads(args.rule.read_text())
    panel=list(csv.DictReader(args.panel.open(encoding="utf-8")))
    expected_hosts={}
    for row in panel:
        ids=[x for x in row["accepted_host_ids"].split(";") if x]
        names=[x for x in row["accepted_host_names"].split(";") if x]
        for i,n in zip(ids,names): expected_hosts[str(i)]=n

    ledgers=[]; occurrence_rows=[]
    for p in sorted(args.shard_dir.rglob("ledger_*.json")):
        ledgers.append(json.loads(p.read_text()))
    for p in sorted(args.shard_dir.rglob("occurrences_*.csv")):
        occurrence_rows.extend(csv.DictReader(p.open(encoding="utf-8")))
    host_meta={}
    for payload in ledgers:
        for row in payload["hosts"]:
            hid=str(row["accepted_host_id"])
            if hid in host_meta: raise RuntimeError(f"duplicate host ledger: {hid}")
            host_meta[hid]=row
    missing=set(expected_hosts)-set(host_meta)
    extra=set(host_meta)-set(expected_hosts)
    if missing or extra:
        raise RuntimeError(f"host ledger mismatch missing={len(missing)} extra={len(extra)}")

    errors=[x for x in host_meta.values() if x["status"]=="REQUEST_ERROR"]
    supported={hid for hid,x in host_meta.items() if x["status"]=="PASS_HOST_OCCURRENCE"}
    completeness=float(rule["insect_host_completeness"]["required_fraction"])
    insect_rows=[]
    for row in panel:
        ids=[x for x in row["accepted_host_ids"].split(";") if x]
        required=int(math.ceil(completeness*len(ids)))
        n=sum(i in supported for i in ids)
        insect_rows.append({
            "species":row["species"],"n_hosts":len(ids),
            "required_supported_hosts":required,
            "supported_hosts":n,
            "occurrence_gate_pass":n>=required,
        })
    retained=[x for x in insect_rows if x["occurrence_gate_pass"]]
    min_insects=int(rule["insect_host_completeness"]["minimum_insects_after_occurrence_gate"])
    if errors:
        status="TECHNICALLY_INCOMPLETE_HISTORICAL_HOST_OCCURRENCE_GATE"
    elif len(retained)>=min_insects:
        status="PASS_TO_PALAEO_HOST_MODEL"
    else:
        status=rule["insect_host_completeness"]["failure_state"]

    # Split only after response-blind occurrence filtering, exactly as prospectively frozen.
    namespace="historical-host-connectivity-panel-v0.1"
    ranked=sorted(
        [x["species"] for x in retained],
        key=lambda s:(hashlib.sha256(f"{namespace}|{s}".encode()).hexdigest(),s),
    )
    ndev=len(ranked)//3
    roles={s:("development" if i<ndev else "confirmatory") for i,s in enumerate(ranked)}
    for x in insect_rows:
        x["role"]=roles.get(x["species"],"")

    args.output_dir.mkdir(parents=True,exist_ok=True)
    occ_path=args.output_dir/"historical_host_occurrences.csv"
    fields=["accepted_host_id","accepted_host_name","source_key","latitude","longitude","priority_rank","priority_sha256"]
    with occ_path.open("w",newline="",encoding="utf-8") as f:
        w=csv.DictWriter(f,fieldnames=fields); w.writeheader()
        w.writerows(sorted(occurrence_rows,key=lambda r:(r["accepted_host_name"],int(r["priority_rank"]))))
    insect_path=args.output_dir/"historical_host_occurrence_insects.csv"
    with insect_path.open("w",newline="",encoding="utf-8") as f:
        w=csv.DictWriter(f,fieldnames=list(insect_rows[0])); w.writeheader(); w.writerows(insect_rows)
    host_path=args.output_dir/"historical_host_occurrence_hosts.csv"
    host_rows=[host_meta[k] for k in sorted(host_meta,key=lambda k:(expected_hosts[k],k))]
    fields_host=sorted({k for r in host_rows for k in r})
    with host_path.open("w",newline="",encoding="utf-8") as f:
        w=csv.DictWriter(f,fieldnames=fields_host); w.writeheader(); w.writerows(host_rows)

    payload={
        "schema":"ttf_genetic_historical_host_occurrence_result_v0.1",
        "status":status,
        "panel_sha256":sha256_path(args.panel),
        "rule_sha256":sha256_path(args.rule),
        "hosts":{
            "total":len(expected_hosts),"supported":len(supported),
            "request_errors":len(errors),
            "failed_or_rejected":len(expected_hosts)-len(supported)-len(errors),
        },
        "insects":{
            "total":len(panel),"passing_completeness":len(retained),
            "minimum_required":min_insects,
            "development":sum(v=="development" for v in roles.values()),
            "confirmatory":sum(v=="confirmatory" for v in roles.values()),
        },
        "outputs_sha256":{
            "occurrences":sha256_path(occ_path),
            "insects":sha256_path(insect_path),
            "hosts":sha256_path(host_path),
        },
        "response_firewall":{
            "fresh_sequence_identity_opened":False,
            "fresh_pairwise_genetic_distances_opened":False,
            "fresh_post_ibd_turnover_computed":False,
            "fresh_beta_host_LGM_computed":False,
        },
    }
    (args.output_dir/"historical_host_occurrence_result.json").write_text(json.dumps(payload,indent=2,sort_keys=True)+"\n")
    print(json.dumps(payload,sort_keys=True))
    return 0 if not errors else 2


if __name__=="__main__":
    raise SystemExit(main())
