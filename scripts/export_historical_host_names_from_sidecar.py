#!/usr/bin/env python3
"""Export canonical WCVP host names from the exact frozen sidecar.

This is transport-only. The scientific host IDs are already frozen.
"""
from __future__ import annotations
import argparse, csv
from pathlib import Path

def split_ids(value: str) -> list[str]:
    return [x.strip() for x in str(value).split(";") if x.strip()]

def main() -> int:
    ap=argparse.ArgumentParser()
    ap.add_argument("--panel",type=Path,required=True)
    ap.add_argument("--native-sidecar",type=Path,required=True)
    ap.add_argument("--output",type=Path,required=True)
    args=ap.parse_args()

    with args.panel.open(newline="",encoding="utf-8") as f:
        panel=list(csv.DictReader(f))
    ids=sorted({hid for row in panel for hid in split_ids(row["accepted_host_ids"])})

    names: dict[str,set[str]]={}
    with args.native_sidecar.open(newline="",encoding="utf-8") as f:
        for row in csv.DictReader(f):
            hid=str(row.get("accepted_plant_name_id","")).strip()
            name=str(row.get("accepted_name","")).strip()
            if hid in ids and name:
                names.setdefault(hid,set()).add(name)

    if set(names)!=set(ids):
        missing=sorted(set(ids)-set(names))
        extra=sorted(set(names)-set(ids))
        raise RuntimeError(f"exact sidecar name coverage mismatch missing={missing} extra={extra}")
    bad={hid:sorted(v) for hid,v in names.items() if len(v)!=1}
    if bad:
        raise RuntimeError(f"exact sidecar canonical-name ambiguity: {bad}")

    args.output.parent.mkdir(parents=True,exist_ok=True)
    with args.output.open("w",newline="",encoding="utf-8") as f:
        w=csv.DictWriter(f,fieldnames=["accepted_host_id","accepted_host_name"])
        w.writeheader()
        for hid in ids:
            w.writerow({"accepted_host_id":hid,"accepted_host_name":next(iter(names[hid]))})
    print(f"hosts={len(ids)}")
    return 0

if __name__=="__main__":
    raise SystemExit(main())
