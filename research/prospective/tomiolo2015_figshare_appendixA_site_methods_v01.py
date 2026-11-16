"""Only original Tomiolo et al 2015 Appendix A site/method text excerpts.

Publisher-linked Figshare collection metadata already qualified. Verify the
exact article DOI, original Appendix A file ID, size, MD5 before parsing.
Extract bounded text *site description and locality context* from appendix A,
explicitly not raw responses, appendix B/C, or fitness effect numbers.
"""
from __future__ import annotations
import hashlib
import html
import json
import re
from html.parser import HTMLParser
from pathlib import Path
from urllib.error import HTTPError,URLError
from urllib.parse import urlsplit
from urllib.request import Request,urlopen

CONTRACT=Path(__file__).with_name("tomiolo2015_figshare_appendixA_site_methods_contract_v01.json")

def canonical_doi(v):
    if not isinstance(v,str):return ""
    return re.sub(r"^(?:https?://(?:dx\.)?doi\.org/|doi:)","",v.strip().lower()).rstrip("/")

def verified_file_url(article:dict,c:dict)->str:
    src=c["source"];authorized=c["authorization"]
    if (article.get("id")!=src["article_id"]
        or canonical_doi(article.get("doi"))!=src["article_doi"]):
        raise ValueError("Not same frozen Figshare original Appendix A article DOI")
    files=article.get("files") or []
    matches=[f for f in files if (isinstance(f,dict)
                and f.get("id")==src["original_file_id"]
                and f.get("name")==src["original_filename"])]
    if len(matches)!=1:raise ValueError("Exact original site Appendix A file absent")
    f=matches[0]
    if (f.get("size")!=src["published_bytes"]
        or str(f.get("computed_md5") or f.get("supplied_md5") or "").lower()!=src["published_md5"]):
        raise ValueError("Appendix A published original bytes/MD5 metadata drift")
    link=f.get("download_url")
    if not isinstance(link,str):raise ValueError("Original Appendix A download URL missing")
    uri=urlsplit(link)
    if (uri.scheme!="https"
        or uri.hostname not in authorized["approved_file_download_hostnames"]
        or not (uri.path.endswith("/files/"+str(src["original_file_id"]))
                or uri.path.endswith("/"+str(src["original_file_id"])))):
        raise ValueError("Appendix A file download is not from official Figshare URL")
    return link

class TextOnly(HTMLParser):
    def __init__(self):
        super().__init__(convert_charrefs=True)
        self.fragments=[]
        self.skip=0
    def handle_starttag(self,name,attrs):
        if name in ("style","script"):self.skip+=1
    def handle_endtag(self,name):
        if name in ("style","script"):self.skip=max(0,self.skip-1)
    def handle_data(self,text):
        if not self.skip and text.strip():self.fragments.append(text)

def source_site_only_excerpts(content:bytes,contract:dict)->dict:
    if (len(content)!=contract["source"]["published_bytes"]
        or hashlib.md5(content).hexdigest()!=contract["source"]["published_md5"]):
        raise ValueError("Unverified original source HTML; no text allowed")
    parser=TextOnly()
    parser.feed(content.decode("utf-8-sig",errors="replace"))
    cleaned=" ".join(parser.fragments)
    cleaned=re.sub(r"\s+"," ",html.unescape(cleaned)).strip()
    allow=contract["authorization"]["allowed_context_words"]
    exclude=contract["authorization"]["exclude_context_if_contains_biological_outcome_terms"]
    pattern=re.compile("|".join(re.escape(k) for k in allow),re.I)
    matches=[]
    seen=set()
    for found in pattern.finditer(cleaned):
        a=max(0,found.start()-70)
        b=min(len(cleaned),found.end()+170)
        excerpt=cleaned[a:b].strip()
        if len(excerpt)>contract["authorization"]["max_context_chars_per_source_match"]:
            excerpt=excerpt[:contract["authorization"]["max_context_chars_per_source_match"]]
        if any(re.search(r"\b"+re.escape(word)+r"\b",excerpt,re.I) for word in exclude):
            continue
        if excerpt.lower() in seen:continue
        seen.add(excerpt.lower())
        matches.append({"site_description_context":excerpt,
                        "original_article_level_experimental_site_context_only":True})
        if len(matches)>=contract["authorization"]["max_relevant_site_context_segments"]:
            break
    return {
      "site_context_excerpts_limited":matches,
      "text_context_count":len(matches),
      "original_html_md5_verified":contract["source"]["published_md5"],
      "original_source_population_points_not_inferred_from_context_alone":True,
      "source_experiment_independence_not_certified":True,
      "numeric_fitness_or_experimental_response_not_extracted":True
    }

def probe(c:dict)->dict:
    auth=c["authorization"]
    if (c["status"]!="FROZEN_BEFORE_ORIGINAL_FIGSHARE_APPENDIX_A_FILE_FETCH"
        or not auth["only_exact_article_and_original_file_metadata_qualified"]
        or not auth["no_appendixB_or_appendixC_fetch"]):
        raise ValueError("Original source-only Appendix A contract drift")
    output={
      "schema":"ttf_Tomiolo2015_Figshare_original_appendixA_source_only_site_methods_v01",
      "original_article_DOI":c["source"]["paper_doi"],
      "original_Figshare_source_DOI":c["source"]["article_doi"],
      "original_appendix_file_id":c["source"]["original_file_id"],
      "original_file_downloaded":False,
      "outcome_values_or_original_data_table_rows_read":False,
      "source_site_population_ID_crosswalk_certified":False,
      "historical_21ka_vs_0BP_climate_effect_fitted":False
    }
    try:
        url=c["source"]["article_metadata_api"]
        if urlsplit(url).hostname!=auth["allowed_article_metadata_hostname"]:
            raise ValueError("Figshare metadata host differs")
        with urlopen(Request(url,headers={"Accept":"application/json",
             "User-Agent":"TTF-original-Tomiolo-appA-only-metadata/2026-11-16"}),timeout=30) as stream:
            meta=stream.read(800001)
        if len(meta)>800000:raise ValueError("Source metadata too large")
        details=json.loads(meta)
        link=verified_file_url(details,c)
        with urlopen(Request(link,headers={"Accept":"*/*",
             "User-Agent":"TTF-original-Tomiolo-AppA-site-methods-only/2026-11-16"}),timeout=45) as response:
            content=response.read(c["source"]["published_bytes"]+1)
        excerpts=source_site_only_excerpts(content,c)
        output.update(excerpts)
        output["original_file_downloaded"]=True
        output["source_file_SHA256"]=hashlib.sha256(content).hexdigest()
        output["status"]="OFFICIAL_APPENDIX_A_SITE_CONTEXT_RETRIEVED_NOT_POPULATION_GPS"
    except (HTTPError,URLError,TimeoutError,ValueError,json.JSONDecodeError,UnicodeError) as e:
        output.update({"status":"OFFICIAL_APPENDIX_A_ACCESS_OR_IDENTITY_NOT_EVALUABLE",
                      "error_type":type(e).__name__,"error_http_code":getattr(e,"code",None)})
    return output

if __name__=="__main__":
    print("TOMIOLO_APPENDIX_A_SITE_ONLY="+json.dumps(
      probe(json.loads(CONTRACT.read_text(encoding="utf-8"))),
      ensure_ascii=False,sort_keys=True))
