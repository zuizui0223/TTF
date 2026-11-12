"""Resolve exact source file receipts behind published Figshare private-share links."""
from __future__ import annotations
import argparse,json,re,urllib.parse,urllib.request,urllib.error
from pathlib import Path
SOURCES={
 "metadata":(25320499,"7a58e7a1b6b3c647b4d2",44757382),
 "linework":(25320052,"1fa446cbb00bad232db3",44755531),
 "code":(25321141,"36bf9c11338c35a9a5e7",44803807)
}
def probe_one(kind,aid,share,first):
  base=f"https://api.figshare.com/v2/articles/{aid}"
  url=base+"?private_link="+share
  req=urllib.request.Request(url,headers={"User-Agent":"TTF-source-receipt/0.1","Accept":"application/json"})
  raw=None;status=None
  try:
    with urllib.request.urlopen(req,timeout=45) as r:
      raw=r.read(1_000_001);status=r.status
      if len(raw)>1_000_000:raise ValueError("metadata oversized")
  except urllib.error.HTTPError as e:status=e.code
  except Exception as e:status="network-error:"+str(e)[:100]
  info={"article_id":aid,"first_file_id_from_page":first,
    "source_share_url":f"https://figshare.com/s/{share}",
    "api_status":status,"api_file_receipts":[],
    "response_geometry_opened":False}
  if status==200 and raw:
    try:
      j=json.loads(raw)
      if int(j.get("id",-1))!=aid:raise ValueError("article id mismatch")
      info["title"]=j.get("title");info["license"]=j.get("license");info["doi"]=j.get("doi")
      for f in j.get("files",[]):
        info["api_file_receipts"].append({k:f.get(k) for k in
           ("id","name","size","supplied_md5","computed_md5","is_link_only","download_url")})
    except ValueError: info["api_parse_error"]="unreadable metadata"
  for file_id in [first]:
    fileurl=f"https://figshare.com/ndownloader/files/{file_id}?private_link={share}"
    try:
      request=urllib.request.Request(fileurl,method="HEAD",headers={"User-Agent":"TTF-source-receipt/0.1"})
      with urllib.request.urlopen(request,timeout=45) as resp:
        info["file_head"]={"status":resp.status,"length":resp.headers.get("Content-Length"),
          "content_type":resp.headers.get("Content-Type"),
          "content_disposition":resp.headers.get("Content-Disposition"),
          "final_domain":urllib.parse.urlparse(resp.url).hostname}
    except urllib.error.HTTPError as e:
      info["file_head"]={"status":e.code}
    except Exception as e:info["file_head"]={"status":"network-error","message":str(e)[:140]}
  return info

def main():
 p=argparse.ArgumentParser();p.add_argument("--out",type=Path,required=True);a=p.parse_args()
 receipts=[probe_one(kind,*items)|{"source_kind":kind} for kind,items in SOURCES.items()]
 result={"schema":"ttf_jensen2024_original_figshare_receipts_v01",
         "status":"FILE_IDENTITIES_ONLY_NO_BREAK_COORDINATES",
         "receipts":receipts,"genetic_break_line_labels_or_coordinates_read":False}
 a.out.parent.mkdir(parents=True,exist_ok=True)
 a.out.write_text(json.dumps(result,indent=2,sort_keys=True)+"\n")
 print(json.dumps([{"kind":x["source_kind"],"api_status":x["api_status"],
   "files":len(x["api_file_receipts"]),"head":x.get("file_head")} for x in receipts]))
if __name__=="__main__":main()
