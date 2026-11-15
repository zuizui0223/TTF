"""Costus OSF nested original-source manifest, no file payloads/biology values.

Roots/folder IDs were observed in public OSF metadata on Actions #43699704360.
Exact original files remain on OSF; never store raw genomic/scent data here.
"""
from __future__ import annotations
import collections
import json
import urllib.error
import urllib.request
from urllib.parse import urljoin

FOLDERS={
  "2025_flower_scent":{
    "node":"2ap4k","folder":"675c524891cac34e4187be6c",
    "name":"Floral Scent Analyses",
  },
  "2025_phylogenies":{
    "node":"2ap4k","folder":"663a31834664da4294ed6b3d",
    "name":"Phylogenies",
  },
  "2022_phylogenomics_pollination":{
    "node":"fkj2x","folder":"62c65b721bb7a500f61f45c4",
    "name":"Genetic mechanisms of convergent pollination syndromes",
  },
}
MAX_PAGES=20
MAX_RECORDS=300
MAX_DEPTH=4

def get_json(url):
    req=urllib.request.Request(url,headers={
       "Accept":"application/vnd.api+json",
       "User-Agent":"TTF-Costus-published-OSF-manifest/1.0",
    })
    with urllib.request.urlopen(req,timeout=35) as resp:
        data=resp.read(2_000_001)
    if len(data)>2_000_000:raise ValueError("OSF metadata oversized")
    return json.loads(data)

def folder_href(record):
    links=(record.get("relationships") or {}).get("files",{}).get("links",{}).get("related")
    if isinstance(links,str):return links
    if isinstance(links,dict):return links.get("href")
    return None

def crawl(folder):
    node=folder["node"];root=folder["folder"]
    queue=collections.deque([(f"https://api.osf.io/v2/nodes/{node}/files/osfstorage/{root}/?page[size]=100",0,folder["name"])])
    seen=set()
    rows=[]
    errors=[]
    while queue and len(seen)<MAX_PAGES and len(rows)<MAX_RECORDS:
        url,depth,prefix=queue.popleft()
        if url in seen:continue
        seen.add(url)
        try:resp=get_json(url)
        except urllib.error.HTTPError as e:
            errors.append({"path":prefix,"http":e.code})
            continue
        items=resp.get("data") or []
        if not isinstance(items,list):raise ValueError("OSF folder response not list")
        for rec in items:
            attrs=rec.get("attributes") or {}
            links=rec.get("links") or {}
            extra=attrs.get("extra") or {}
            name=attrs.get("name")
            if not isinstance(name,str) or not name:raise ValueError("OSF file entry without name")
            kind=attrs.get("kind")
            obj={
                "path":prefix+"/"+name,
                "name":name,
                "kind":kind,
                "size":attrs.get("size"),
                "id":rec.get("id"),
                "sha256":extra.get("hashes",{}).get("sha256") if isinstance(extra,dict) else None,
                "md5":extra.get("hashes",{}).get("md5") if isinstance(extra,dict) else None,
                "download":links.get("download") if kind=="file" else None,
            }
            rows.append(obj)
            if kind=="folder" and depth<MAX_DEPTH:
                sub=folder_href(rec)
                if sub:
                    queue.append((urljoin("https://api.osf.io",sub),depth+1,prefix+"/"+name))
        next_page=(resp.get("links") or {}).get("next")
        if next_page:
            queue.appendleft((urljoin("https://api.osf.io",next_page),depth,prefix))
    if queue: errors.append({"stop":"PAGE_OR_RECORD_CAP","pages":len(seen),"entries":len(rows)})
    return {"status":"OSF_NESTED_FILE_METADATA_ONLY","root":folder,"files":rows,
            "listed_files":sum(x["kind"]=="file" for x in rows),
            "listed_folders":sum(x["kind"]=="folder" for x in rows),
            "metadata_errors":errors,
            "no_data_file_downloaded":True}

def main():
    results={name:crawl(group) for name,group in FOLDERS.items()}
    print(json.dumps({"status":"COSTUS_MULTIPAPER_OSF_SOURCE_INVENTORY",
        "sources":results,"outcome_values_read":False,
        "not_a_species_matched_analysis":True},indent=2,sort_keys=True,ensure_ascii=False))

if __name__=="__main__":main()
