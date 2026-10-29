"""Response-blind file listing for independent GenDivRange biological validation."""
from __future__ import annotations
import argparse
import json
import re
import urllib.request
from pathlib import Path

ARTICLE_ID=28407989
API_URL=f"https://api.figshare.com/v2/articles/{ARTICLE_ID}"
VALID_EXT={".csv",".tsv",".xlsx",".xls",".zip",".rds",".rdata",".parquet"}

def normalize_metadata(raw:dict)->dict:
    if not isinstance(raw,dict) or int(raw.get("id",-1))!=ARTICLE_ID:
        raise ValueError("unexpected Figshare article identity")
    files=raw.get("files")
    if not isinstance(files,list):raise ValueError("missing files list")
    records=[]
    for f in files:
        if not isinstance(f,dict) or any(k not in f for k in ("id","name","size","download_url")):
            raise ValueError("unexpected file record")
        name=str(f["name"])
        ext=Path(name).suffix.lower()
        records.append({"id":int(f["id"]),"filename":re.sub(r"[^a-zA-Z0-9_.-]+","_",name)[:160],
            "size_bytes":int(f["size"]),"ext":ext,"data_candidate":ext in VALID_EXT,
            "md5_expected":f.get("supplied_md5"),"download_url":str(f["download_url"])})
    return {"schema":"gendivrange_metadata_census_v01",
        "status":"FILE_LIST_ONLY_NO_GENETIC_RESPONSE",
        "figshare_article_id":ARTICLE_ID,"article_doi":str(raw.get("doi","")),
        "file_count":len(records),"files":records,"biology_response_values_accessed":False}

def fetch_metadata()->dict:
    req=urllib.request.Request(API_URL,headers={"User-Agent":"TTF-independent-response-blind/0.1"})
    with urllib.request.urlopen(req,timeout=40) as response:
        data=response.read(2_000_001)
    if len(data)>2_000_000:raise ValueError("unusually large metadata")
    return json.loads(data)

def main()->None:
    p=argparse.ArgumentParser()
    p.add_argument("--output",type=Path,required=True)
    p.add_argument("--input-metadata-json",type=Path)
    args=p.parse_args()
    raw=json.loads(args.input_metadata_json.read_text()) if args.input_metadata_json else fetch_metadata()
    result=normalize_metadata(raw)
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(result,indent=2,sort_keys=True)+"\n")
    print(json.dumps({"status":result["status"],"file_count":result["file_count"],
        "file_names":[x["filename"] for x in result["files"]]}))
if __name__=="__main__":main()
