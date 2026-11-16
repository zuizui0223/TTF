"""Official Tomiolo 2015 source-only Figshare collection/file manifest.

Original paper: Ecology DOI 10.1890/14-1445.1.
Original publisher-linked Figshare collection DOI
10.6084/m9.figshare.c.3307653.v1.
Only official JSON records and attachment *filenames, sizes and digests* are
retrieved. No files, tables, experimental fitness or coordinates are read.
"""
from __future__ import annotations
import json
import re
from pathlib import Path
from urllib.parse import urlsplit
from urllib.request import Request,urlopen
from urllib.error import HTTPError,URLError

CONTRACT=Path(__file__).with_name("tomiolo2015_source_geo_figshare_manifest_contract_v01.json")

def get_official_metadata(url:str,policy:dict):
    parsed=urlsplit(url)
    if parsed.scheme!="https" or parsed.hostname not in policy["only_public_official_api_hosts"]:
        raise ValueError("Non-official Figshare API URL prohibited")
    req=Request(url,headers={
        "Accept":"application/json",
        "User-Agent":"TTF-Tomiolo2015-source-only-manifest/2026-11-16"})
    with urlopen(req,timeout=30) as response:
        data=response.read(policy["max_metadata_response_bytes"]+1)
    if len(data)>policy["max_metadata_response_bytes"]:
        raise ValueError("Official metadata size exceeds frozen cap")
    return json.loads(data)

def doi_norm(value):
    if not isinstance(value,str):return ""
    return re.sub(r"^(?:https?://(?:dx\.)?doi\.org/|doi:)","",value.strip().lower()).rstrip("/")

def valid_collection(collection:dict,contract:dict):
    target=contract["target"]
    if str(collection.get("id"))!=str(contract["figshare_source_manifest_policy"]["exact_collection_id_required"]):
        raise ValueError("Wrong source collection id")
    version_doi=doi_norm(collection.get("doi"))
    source_doi=doi_norm(target["figshare_collection_doi"])
    if version_doi!=source_doi and version_doi!=source_doi.removesuffix(".v1"):
        raise ValueError("Figshare collection original DOI mismatch")
    return version_doi

def files_from_article(article:dict,policy:dict)->list[dict]:
    files=article.get("files") or []
    if not isinstance(files,list) or len(files)>policy["max_files_per_article"]:
        raise ValueError("Original public Figshare file listing outside allowed bound")
    output=[]
    for f in files:
        if not isinstance(f,dict):
            raise ValueError("Unexpected Figshare file object")
        name=f.get("name")
        if not isinstance(name,str) or len(name)>255:
            raise ValueError("Figshare original filename invalid")
        output.append({
            "original_file_id":f.get("id"),
            "original_file_name":name,
            "published_size_bytes":f.get("size"),
            "original_content_checksum_from_metadata":f.get("computed_md5") or f.get("supplied_md5"),
            "no_source_file_body_opened":True
        })
    return sorted(output,key=lambda x:x["original_file_name"])

def run()->dict:
    c=json.loads(CONTRACT.read_text(encoding="utf-8"))
    policy=c["figshare_source_manifest_policy"]
    if c["status"]!="FROZEN_BEFORE_TOMIOLO_FIGSHARE_PUBLIC_DATA_MANIFEST_QUERY":
        raise ValueError("No frozen source-only Tomiolo manifest design")
    if policy["data_files_downloaded"] or policy["max_article_detail_requests"]!=12:
        raise ValueError("External file access contract drifted")
    attempts=[]
    result={
        "schema":"ttf_Tomiolo2015_original_Figshare_collection_source_only_file_manifest_v01",
        "original_paper_doi":c["target"]["original_paper_doi"],
        "original_data_collection_doi":c["target"]["figshare_collection_doi"],
        "target_original_study_label":c["target"]["hargreaves_original_label"],
        "Hargreaves_six_source_ID_keys_map_to_only_two_distinct_garden_proxy_points":True,
        "external_thesis_sites_M_and_SA_with_whole_minute_coordinates":True,
        "source_experiment_primary_independence_certified":False,
        "original_per_source_population_GPS_certified":False,
        "all_data_file_bodies_or_biological_outcome_values_opened":False,
        "figshare_public_collection_metadata":None,
        "figshare_public_article_manifest":[],
        "attempts":attempts
    }
    try:
        col=get_official_metadata(c["target"]["figshare_collection_api"],policy)
        doi=valid_collection(col,c)
        result["figshare_public_collection_metadata"]={
            "id":col.get("id"),"doi":doi,"title":str(col.get("title") or "")[:200],
            "public_url":col.get("url_public_html") or col.get("url_public_api")
        }
        attempts.append({"phase":"exact_source_collection",
                         "status":"DOI_AND_ID_EXACT_VALIDATED"})
    except (HTTPError,URLError,TimeoutError,ValueError,json.JSONDecodeError) as exc:
        attempts.append({"phase":"exact_source_collection","status":"ACCESS_OR_IDENTITY_UNQUALIFIED",
                         "error_type":type(exc).__name__,"http":getattr(exc,"code",None)})
        result["status"]="ORIGINAL_FIGSHARE_COLLECTION_METADATA_NOT_EVALUABLE"
        return result
    try:
        lst=get_official_metadata(c["target"]["figshare_collection_articles_api"],policy)
        if not isinstance(lst,list) or len(lst)>policy["max_collection_articles"]:
            raise ValueError("Original collection article listing outside frozen bounds")
        attempts.append({"phase":"public_collection_articles","status":"LIST_RETURNED","count":len(lst)})
        partial=len(lst)>policy["max_article_detail_requests"]
        accepted=lst[:policy["max_article_detail_requests"]]
        for v in accepted:
            articleid=v.get("id") if isinstance(v,dict) else None
            if not isinstance(articleid,int) or articleid<1:
                raise ValueError("Unexpected Figshare article metadata id")
            try:
                details=get_official_metadata(
                    f"https://api.figshare.com/v2/articles/{articleid}",policy)
                if details.get("id")!=articleid:
                    raise ValueError("Figshare original article identifier mismatch")
                result["figshare_public_article_manifest"].append({
                    "id":articleid,
                    "title":str(details.get("title") or "")[:160],
                    "article_doi":details.get("doi"),
                    "url_public_html":details.get("url_public_html"),
                    "original_file_metadata":files_from_article(details,policy),
                    "response_values_or_original_data_rows_read":False
                })
            except (HTTPError,URLError,TimeoutError,ValueError,json.JSONDecodeError) as exc:
                result["figshare_public_article_manifest"].append({
                    "id":articleid,"status":"ARTICLE_METADATA_NOT_EVALUABLE",
                    "error_type":type(exc).__name__,"http":getattr(exc,"code",None)
                })
        result["status"]="PUBLIC_FIGSHARE_SOURCE_MANIFEST_PARTIAL_NO_DATA_FILES" if partial else "PUBLIC_FIGSHARE_SOURCE_MANIFEST_LISTED_NO_DATA_FILES"
        result["articles_in_public_collection_listing"]=len(lst)
        result["articles_inspected_up_to_precommitted_limit"]=len(accepted)
        result["metadata_discovery_complete"]=not partial
    except (HTTPError,URLError,TimeoutError,ValueError,json.JSONDecodeError) as exc:
        attempts.append({"phase":"collection_articles","status":"ACCESS_OR_METADATA_ERROR",
                         "error_type":type(exc).__name__,"http":getattr(exc,"code",None)})
        result["status"]="PUBLIC_FIGSHARE_SOURCE_COLLECTION_IDENTIFIED_FILES_NOT_EVALUABLE"
    return result

if __name__=="__main__":
    print("TOMIOLO2015_FIGSHARE_SOURCE_MANIFEST="+json.dumps(run(),ensure_ascii=False,sort_keys=True))
