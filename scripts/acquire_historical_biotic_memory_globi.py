#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
from pathlib import Path
import time
from urllib.parse import urlencode
from urllib.request import Request, urlopen


EXPECTED_RULE_SCHEMA="ttf_historical_biotic_memory_globi_capture_rule_v0.1"


def sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def sha256_path(path: Path) -> str:
    h=hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda:handle.read(1<<20),b""):
            h.update(chunk)
    return h.hexdigest()


def read_species(path: Path) -> list[str]:
    rows=[line.strip() for line in path.read_text().splitlines() if line.strip()]
    if len(rows)!=len(set(rows)):
        raise RuntimeError("duplicate species in frozen panel")
    return rows


def fetch_bytes(url: str, *, retries: int, timeout: int) -> bytes:
    error: Exception | None=None
    for attempt in range(int(retries)):
        try:
            req=Request(
                url,
                headers={
                    "User-Agent":"TTF-historical-biotic-memory-v0.1/response-blind-capture",
                    "Accept":"text/csv,*/*;q=0.8",
                },
            )
            with urlopen(req,timeout=int(timeout)) as response:
                body=response.read()
                if int(getattr(response,"status",200))!=200:
                    raise RuntimeError(f"HTTP status {response.status}")
                return body
        except Exception as exc:
            error=exc
            if attempt+1<int(retries):
                time.sleep(2**attempt)
    raise RuntimeError(f"GloBI request failed after {retries} attempts: {error}")


def csv_shape(body: bytes) -> tuple[list[str],int]:
    text=body.decode("utf-8-sig")
    reader=csv.reader(io.StringIO(text))
    try:
        header=next(reader)
    except StopIteration:
        return [],0
    return [str(x) for x in header],sum(1 for _ in reader)


def main() -> int:
    ap=argparse.ArgumentParser()
    ap.add_argument("--rule",type=Path,required=True)
    ap.add_argument("--species",type=Path,required=True)
    ap.add_argument("--start",type=int,required=True)
    ap.add_argument("--stop",type=int,required=True)
    ap.add_argument("--output-dir",type=Path,required=True)
    args=ap.parse_args()

    rule=json.loads(args.rule.read_text())
    if rule.get("schema")!=EXPECTED_RULE_SCHEMA:
        raise RuntimeError("unexpected GloBI capture rule")
    if rule.get("status")!="FROZEN_BEFORE_GLOBI_RAW_CAPTURE":
        raise RuntimeError("GloBI capture rule is not frozen")
    if any(bool(v) for v in rule["response_firewall"].values()):
        raise RuntimeError("genetic response firewall is open")

    expected=rule["species"]
    if sha256_path(args.species)!=expected["sha256"]:
        raise RuntimeError("frozen species file SHA drift")
    species=read_species(args.species)
    if len(species)!=int(expected["count"]):
        raise RuntimeError("frozen species count drift")
    if not (0<=args.start<args.stop<=len(species)):
        raise RuntimeError("invalid shard bounds")

    out=args.output_dir
    raw=out/"raw"
    raw.mkdir(parents=True,exist_ok=True)
    entries=[]
    headers_union=set()
    base=rule["source"]["api_base"]
    retries=int(rule["capture"]["retries_per_request"])
    timeout=int(rule["capture"]["timeout_seconds"])

    for index in range(args.start,args.stop):
        name=species[index]
        for direction in rule["capture"]["directions"]:
            url=base+"?"+urlencode({
                "includeObservations":rule["source"]["query_parameter_includeObservations"],
                direction:name,
            })
            body=fetch_bytes(url,retries=retries,timeout=timeout)
            header,nrows=csv_shape(body)
            headers_union.update(header)
            slug=f"{index:03d}_{direction}"
            target=raw/f"{slug}.csv"
            target.write_bytes(body)
            entries.append({
                "species_index":index,
                "species":name,
                "direction":direction,
                "query_url":url,
                "raw_file":str(target.relative_to(out)),
                "raw_sha256":sha256_bytes(body),
                "bytes":len(body),
                "rows":int(nrows),
                "header":header,
            })

    payload={
        "schema":"ttf_historical_biotic_memory_globi_raw_capture_shard_v0.1",
        "status":"RAW_GLOBI_CAPTURE_COMPLETE",
        "bounds":[args.start,args.stop],
        "species_count":args.stop-args.start,
        "queries":len(entries),
        "entries":entries,
        "headers_union":sorted(headers_union),
        "response_firewall":rule["response_firewall"],
    }
    out.mkdir(parents=True,exist_ok=True)
    manifest=out/f"globi_capture_{args.start:03d}_{args.stop:03d}.json"
    manifest.write_text(json.dumps(payload,indent=2,sort_keys=True)+"\n")
    print(json.dumps({
        "status":payload["status"],
        "bounds":payload["bounds"],
        "queries":payload["queries"],
        "rows_total":sum(x["rows"] for x in entries),
        "nonempty_queries":sum(x["rows"]>0 for x in entries),
        "manifest_sha256":sha256_path(manifest),
    },sort_keys=True))
    return 0


if __name__=="__main__":
    raise SystemExit(main())
