"""Costus OSF 2025 scent and 2022 phylogenomic original file manifests.

Query public OSF nodes/file metadata only. Do not retrieve raw scent,
transcriptomes, genomic reads, pollinator values, or experimental outcomes.
"""
from __future__ import annotations
import json
import urllib.error
import urllib.request
from urllib.parse import urljoin

NODES={
  "scent_2025":"2ap4k",
  "phylogenomics_2022":"fkj2x",
}
BASE="https://api.osf.io/v2/"
MAX_ENTRIES=200
MAX_PAGES=10


def get_json(url:str)->dict:
    req=urllib.request.Request(url,headers={
        "Accept":"application/vnd.api+json",
        "User-Agent":"TTF-CCOC-Costus-public-metadata/1.0",
    })
    with urllib.request.urlopen(req,timeout=35) as response:
        data=response.read(2_000_001)
    if len(data)>2_000_000:
        raise ValueError("OSF metadata response over cap")
    return json.loads(data)


def meta_file(record:dict)->dict:
    attrs=record.get("attributes") or {}
    links=record.get("links") or {}
    extra=attrs.get("extra") or {}
    return {
        "id":record.get("id"),
        "name":attrs.get("name"),
        "kind":attrs.get("kind"),
        "size":attrs.get("size"),
        "provider":attrs.get("provider"),
        "extra_hashes":extra.get("hashes") if isinstance(extra,dict) else None,
        "path":attrs.get("path"),
        "download_link":links.get("download"),
        "related_children":(record.get("relationships") or {}).get("files",{}).get("links",{}).get("related",{}).get("href"),
    }


def node_inventory(node_id:str)->dict:
    info=get_json(f"{BASE}nodes/{node_id}/")
    a=info.get("data") or {}
    if a.get("id")!=node_id:
        raise ValueError("OSF node identity mismatch")
    attrs=a.get("attributes") or {}
    index={
        "node_id":node_id,
        "title":attrs.get("title"),
        "public":attrs.get("public"),
        "url":f"https://osf.io/{node_id}/",
        "status":"METADATA_ONLY",
    }
    files_url=f"{BASE}nodes/{node_id}/files/"
    providers=get_json(files_url).get("data") or []
    provider_names=[p.get("attributes",{}).get("name") or p.get("id") for p in providers]
    rows=[]
    errors=[]
    current=f"{BASE}nodes/{node_id}/files/osfstorage/?page[size]=100"
    seen=set()
    for _ in range(MAX_PAGES):
        if not current or current in seen:break
        seen.add(current)
        try:
            resp=get_json(current)
        except urllib.error.HTTPError as e:
            errors.append({"url":current,"http":e.code})
            break
        items=resp.get("data") or []
        if not isinstance(items,list):
            raise ValueError("OSF storage payload not array")
        rows += [meta_file(x) for x in items]
        if len(rows)>=MAX_ENTRIES:
            rows=rows[:MAX_ENTRIES]
            errors.append({"stop":"MAX_ENTRIES"})
            break
        current=(resp.get("links") or {}).get("next")
        if current and not current.startswith("https://"):
            current=urljoin(BASE,current)
    index.update({
        "storage_providers":provider_names,
        "osfstorage_top_level_metadata":rows,
        "metadata_errors":errors,
        "top_level_entry_count":len(rows),
        "no_file_payload_downloaded":True,
    })
    return index


def main():
    result={}
    for label,node in NODES.items():
        try:
            result[label]=node_inventory(node)
        except urllib.error.HTTPError as e:
            result[label]={"node_id":node,"status":"PUBLIC_OSF_METADATA_HTTP_BLOCKED",
                           "http":e.code,"no_file_payload_downloaded":True}
    print(json.dumps({"status":"COSTUS_OSF_PUBLIC_METADATA_ONLY",
                      "nodes":result,"no_biological_outcomes_read":True,
                      "none_of_the_listed_files_have_been_downloaded":True},
                     indent=2,sort_keys=True,ensure_ascii=False))


if __name__=="__main__":main()
