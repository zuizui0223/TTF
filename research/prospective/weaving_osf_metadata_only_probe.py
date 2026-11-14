"""Discover published Weaving 2022 OSF file metadata only, never response rows.

This source-only probe lists file IDs, reported checksums, names and sizes from
the public OSF API. It does NOT download biological effect-size datasets, fit
plasticity models, sample climate, read genetic outcomes or alter TTF freezes.
"""
from __future__ import annotations
import argparse
import json
import urllib.request
from urllib.parse import urlparse
from pathlib import Path

HERE=Path(__file__).resolve().parent
PLAN=HERE/"weaving_2022_source_metadata_plan_v01.json"


def safe_osf_url(url):
    parts=urlparse(url)
    if parts.scheme!="https" or parts.netloc!="api.osf.io" or not parts.path.startswith("/v2/"):
        raise ValueError("OSF API metadata only; unsupported URL")
    return url


def metadata_get(url):
    req=urllib.request.Request(safe_osf_url(url),headers={
        "Accept":"application/vnd.api+json",
        "User-Agent":"TTF-phenotype-source-metadata-only-v01",
    })
    with urllib.request.urlopen(req,timeout=40) as r:
        return json.load(r)


def extract_entry(obj):
    a=obj.get("attributes") or {}
    links=obj.get("links") or {}
    rel=obj.get("relationships") or {}
    # DO NOT open links.download or read response content.
    out={
        "id":str(obj.get("id") or ""),
        "name":str(a.get("name") or ""),
        "kind":str(a.get("kind") or ""),
        "size_bytes":a.get("size"),
        "date_modified":a.get("date_modified"),
        "materialized_path":a.get("materialized_path"),
        "checksum_metadata":a.get("extra") if isinstance(a.get("extra"),dict) else None,
    }
    nxt=None
    if out["kind"]=="folder":
        files=rel.get("files") or {}
        more=files.get("links") or {}
        rr=more.get("related") or {}
        if isinstance(rr,dict):
            nxt=rr.get("href")
        elif isinstance(rr,str):
            nxt=rr
        if not nxt:
            api=links.get("self")
            if api and safe_osf_url(api):
                nxt=api.rstrip("/")+"/children/"
    return out,nxt


def crawl(root,fetch=metadata_get):
    seen=set()
    queues=[(root,0)]
    files=[]
    folders=[]
    max_depth=4
    max_pages=50
    while queues:
        link,depth=queues.pop(0)
        if not link or link in seen:
            continue
        safe_osf_url(link)
        if depth>max_depth or len(seen)>=max_pages:
            raise ValueError("Public OSF source listing exceeds frozen safety limits")
        seen.add(link)
        data=fetch(link)
        for obj in data.get("data") or []:
            item,next_folder=extract_entry(obj)
            if item["kind"]=="folder":
                folders.append(item)
                if next_folder:
                    queues.append((next_folder,depth+1))
            elif item["kind"]=="file":
                files.append(item)
            else:
                folders.append(item)
        nextpage=(data.get("links") or {}).get("next")
        if nextpage:
            queues.append((nextpage,depth))
    files.sort(key=lambda z:(z["materialized_path"] or "",z["name"]))
    return {"files":files,"folders":folders,"metadata_pages_retrieved":len(seen)}


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--output",type=Path,required=True)
    args=p.parse_args()
    plan=json.loads(PLAN.read_text())
    if plan["status"]!="PRE_BIOLOGICAL_RESPONSE_METADATA_DISCOVERY_ONLY":
        raise ValueError("Unfrozen metadata-only OSF plan")
    root=plan["source"]["osf_files_api_base"]
    report={"schema":"ttf_weaving_2022_public_osf_source_metadata_only_v01",
        "target_publication_species":plan["source"]["publication_species"],
        "source_project":plan["source"]["osf_public_url"],
        "osf_metadata_api":root,
        "effect_size_rows_downloaded_or_inspected":False,
        "source_coordinates_opened":False,
        "history_climate_extracted":False,
        "genetic_outcomes_opened":False,
        "biological_predictive_test_performed":False}
    try:
        result=crawl(root)
        report.update(result)
        report["source_file_count"]=len(result["files"])
        report["status"]="OSF_FILE_METADATA_DISCOVERED_NO_BIOLOGICAL_RESPONSE_OPENED"
    except Exception as e:
        # A source-access failure does not become an ecological null.
        report["status"]="SOURCE_METADATA_UNAVAILABLE_NOT_A_BIOLOGICAL_RESULT"
        report["error_type"]=type(e).__name__
        report["error_message"]=str(e)[:400]
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(report,sort_keys=True,indent=2)+"\n")
    print(json.dumps(report,sort_keys=True))


if __name__=="__main__":
    main()
