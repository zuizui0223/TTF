"""Metadata-only public source discovery. No genetic break coordinates or labels."""
import argparse,json,re,urllib.parse,urllib.request,urllib.error
from html.parser import HTMLParser
from pathlib import Path

LINKS={
  "metadata":"https://figshare.com/s/7a58e7a1b6b3c647b4d2",
  "linework":"https://figshare.com/s/1fa446cbb00bad232db3",
  "code":"https://figshare.com/s/36bf9c11338c35a9a5e7"
}
API="https://api.figshare.com/v2/articles/search"
class Meta(HTMLParser):
  def __init__(self):super().__init__();self.canonical=[];self.title=[];self.description=[]
  def handle_starttag(self,tag,attrs):
    v=dict(attrs)
    if tag=="link" and v.get("rel")=="canonical":self.canonical.append(v.get("href",""))
    if tag=="meta" and v.get("property") in ("og:title","og:description"):
      (self.title if v["property"]=="og:title" else self.description).append(v.get("content",""))

def safe_request(url, payload=None, limit=1_000_000):
  parts=urllib.parse.urlparse(url)
  if parts.scheme!="https" or parts.hostname not in ("figshare.com","api.figshare.com"):
    raise ValueError("unapproved source hostname")
  headers={"User-Agent":"TTF-mammal-break-source-preflight-v01",
           "Accept":"application/json" if payload is not None else "text/html,application/json"}
  body=json.dumps(payload).encode() if payload is not None else None
  if body is not None:headers["Content-Type"]="application/json"
  req=urllib.request.Request(url,data=body,headers=headers)
  try:
    with urllib.request.urlopen(req,timeout=35) as r:
      data=r.read(limit+1)
      if len(data)>limit:raise ValueError("source returned oversized metadata")
      return {"status":r.status,"final_url":r.url,"raw":data}
  except urllib.error.HTTPError as e:
    return {"status":e.code,"final_url":e.url,"raw":e.read(2048)}
  except Exception as e:
    return {"status":"NETWORK_ERROR","final_url":url,"error":str(e)[:200],"raw":b""}

def search_articles(term):
  response=safe_request(API,{"search_for":term,"page_size":60,"order":"published_date","order_direction":"desc"})
  out={"query":term,"status":response["status"]}
  try:
    j=json.loads(response["raw"]);items=j if isinstance(j,list) else j.get("items",[])
    out["hits"]=[{k:a.get(k) for k in ("id","title","doi","published_date","url_public_html")}
                 for a in items[:60] if isinstance(a,dict)]
  except (ValueError,TypeError,AttributeError):
    out["hits"]=[]
    out["error"]="No searchable JSON results"
  return out
def probe_share(name,url):
  r=safe_request(url)
  meta=Meta()
  if isinstance(r["status"],int) and r["status"]==200:
    meta.feed(r["raw"].decode("utf-8","replace"))
  p=urllib.parse.urlparse(r["final_url"])
  match=re.findall(r"/(?:articles|collections|projects)/[^'\"]+/(\d+)",p.path)
  return {"source_type":name,"initial_url":url,
          "http_status":r["status"],"final_url":r["final_url"],
          "page_title":meta.title[:1],"canonical":meta.canonical[:1],
          "description":meta.description[:1],"article_id_candidates":match,
          "network_error":r.get("error"),
          "response_bytes":len(r["raw"]),"genetic_response_opened":False}
def main():
  a=argparse.ArgumentParser();a.add_argument("--out",type=Path,required=True)
  args=a.parse_args()
  results={"schema":"ttf_jensen2024_genetic_break_source_census_v01",
           "status":"METADATA_DISCOVERY_ONLY","share_probes":[probe_share(k,v) for k,v in LINKS.items()],
           "article_searches":[search_articles(s) for s in
                  ("Geographic barriers but not life history traits shape the phylogeography",
                  "Alex Jensen North American mammal phylogeography",
                  "phylogeographic breaks Jensen 2024")],
           "numeric_or_spatial_genetic_break_labels_opened":False}
  args.out.parent.mkdir(parents=True,exist_ok=True)
  args.out.write_text(json.dumps(results,indent=2,sort_keys=True)+"\n")
  print(json.dumps({"status":results["status"],
      "share_probes":[{k:v for k,v in x.items() if k in ("source_type","http_status","final_url","article_id_candidates")}
                     for x in results["share_probes"]],
      "search_hit_counts":[len(x.get("hits",[])) for x in results["article_searches"]]}))
if __name__=="__main__":main()
