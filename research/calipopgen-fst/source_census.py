"""Discover CaliPopGen Dataset 2 source files, metadata only; never inspect Fst."""
import argparse
import json
import urllib.request
from pathlib import Path
API="https://api.figshare.com/v2"
COLLECTION=5518389
def get_json(url):
    if not url.startswith(API+"/"): raise ValueError("URL must use Figshare API")
    request=urllib.request.Request(url,headers={"User-Agent":"TTF-CaliPopGen-FST-source-census/0.1"})
    with urllib.request.urlopen(request,timeout=45) as response:
        content=response.read(2_000_001)
    if len(content)>2_000_000:raise ValueError("metadata size exceeds cap")
    return json.loads(content)
def choose_article(items):
    if isinstance(items,dict):items=items.get("items",items.get("articles",[]))
    hits=[a for a in items if "calipopgen" in str(a.get("title","")).lower()
          and "dataset 2" in str(a.get("title","")).lower()
          and "pairwise" in str(a.get("title","")).lower()]
    if len(hits)!=1:raise ValueError(f"Expected one distinct Dataset 2, found {len(hits)}")
    return int(hits[0]["id"])
def normalize_detail(detail):
    title=str(detail.get("title",""))
    if "calipopgen" not in title.lower() or "pairwise" not in title.lower():
        raise ValueError("unexpected Figshare article")
    files=[]
    for f in detail.get("files",[]):
        files.append({"id":int(f["id"]),"name":str(f["name"]),
                      "size_bytes":int(f["size"]),"md5":f.get("supplied_md5") or f.get("computed_md5"),
                      "is_link_only":bool(f.get("is_link_only",False))})
    if not files:raise ValueError("no source files")
    return {"schema":"ttf_calipopgen_fst_figshare_metadata_v01",
            "status":"METADATA_ONLY_NO_FST_VALUES",
            "article_id":int(detail["id"]),
            "article_title":title,
            "doi":str(detail.get("doi","")),
            "version":detail.get("version"),"license":detail.get("license"),
            "files":files,
            "Fst_numeric_opened":False}
def main():
    a=argparse.ArgumentParser()
    a.add_argument("--output",required=True,type=Path)
    args=a.parse_args()
    listing=get_json(f"{API}/collections/{COLLECTION}/articles")
    id=choose_article(listing)
    data=normalize_detail(get_json(f"{API}/articles/{id}"))
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(data,indent=2,sort_keys=True)+"\n")
    print(json.dumps({"status":data["status"],"article_id":id,"files":
                      [{"name":v["name"],"size_bytes":v["size_bytes"]} for v in data["files"]]}))
if __name__=="__main__":main()
