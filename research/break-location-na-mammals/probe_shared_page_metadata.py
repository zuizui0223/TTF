"""Extract publication Figshare link target IDs and file receipts; no biology geometry."""
from __future__ import annotations
import argparse,collections,html,json,re,urllib.request
from html.parser import HTMLParser
from pathlib import Path
from source_probe import LINKS

class Sniff(HTMLParser):
 def __init__(self):super().__init__();self.in_script=False;self.script_id=None;self.scripts=[];self.links=[]
 def handle_starttag(self,tag,attrs):
  a=dict(attrs)
  if tag=="script":
   self.in_script=True;self.script_id=(a.get("id"),a.get("type"))
  if tag in ("a","link") and "href" in a:
   h=a["href"]
   if re.search(r"(figshare|ndownloader|/files/|/articles/|\.zip|\.shp|\.gpkg|\.csv|\.rds)",h,re.I):
    self.links.append(h)
 def handle_endtag(self,tag):
  if tag=="script":self.in_script=False;self.script_id=None
 def handle_data(self,data):
  if self.in_script and len(data)>25:self.scripts.append((self.script_id,data))
def fetch(url):
 req=urllib.request.Request(url,headers={"User-Agent":"TTF-Jensen2024-sourcelist-v01"})
 with urllib.request.urlopen(req,timeout=50) as resp:return resp.read(800000).decode("utf-8","replace")
def scan(token,body):
 sniff=Sniff();sniff.feed(body)
 samples={}
 pats={
 "article_id":r'(?:articleId|article_id|itemId|item_id|articleID)["\\\s:=]+(\d{5,12})',
 "article_link":r'(?:articles|filesets|datasets)/(?:[^\s"\\<>]*?/)?(\d{6,12})',
 "figshare_file_url":r'(?:ndownloader\.figshare\.com/files/)(\d+)',
 "file_id":r'(?:fileId|file_id)["\\\s:=]+(\d{5,12})',
 "figshare_data_api":r'api\.figshare\.com/v\d/[^"\s<>]+',
 "figshare_download":r'(?:download_url|downloadUrl)["\\\s:]+.{0,130}'
 }
 for key,pat in pats.items():
  matches=re.findall(pat,body,re.I)
  samples[key]=list(dict.fromkeys(matches))[:25]
 tags=[]
 for meta,src in sniff.scripts:
  sid,typ=meta or (None,None)
  if sid in ("__NEXT_DATA__","__NUXT_DATA__") or typ=="application/ld+json":
   try:
    obj=json.loads(html.unescape(src))
    txt=json.dumps(obj)[:50000]
    ids=list(dict.fromkeys(re.findall(r'(?:"(?:article_id|articleId|id)"\s*:\s*)(\d{6,12})',txt)))
    tags.append({"script_id":sid,"type":typ,"top_keys":list(obj.keys())[:20] if isinstance(obj,dict) else None,
      "object_identifiers":ids[:25],"raw_bytes":len(src)})
   except Exception:tags.append({"script_id":sid,"type":typ,"raw_bytes":len(src),"parse":False})
 return {"name":token,"source_page_bytes":len(body),"embedded_script_count":len(sniff.scripts),
  "matched_metadata_patterns":samples,
  "structured_scripts":tags,
  "candidate_hrefs":list(dict.fromkeys(sniff.links))[:40],
  "geographic_break_label_coordinates_read":False}
def main():
 p=argparse.ArgumentParser();p.add_argument("--output",type=Path,required=True);a=p.parse_args()
 res={"schema":"ttf_jensen_figshare_shared_page_metadata_v01",
  "status":"SOURCE_PAGE_DIAGNOSTIC_ONLY",
  "shares":[scan(k,fetch(url)) for k,url in LINKS.items()],
  "genetic_break_coords_not_opened":True}
 a.output.parent.mkdir(parents=True,exist_ok=True)
 a.output.write_text(json.dumps(res,indent=2,sort_keys=True)+"\n")
 print(json.dumps([{"name":x["name"],"script_count":x["embedded_script_count"],
   "pattern_counts":{k:len(v) for k,v in x["matched_metadata_patterns"].items()},
   "structured_scripts":x["structured_scripts"][:2]} for x in res["shares"]]))
if __name__=="__main__":main()
