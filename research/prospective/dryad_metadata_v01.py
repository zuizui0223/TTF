"""Public Dryad v2 catalog for new transplant studies. Does NOT read fitness rows."""
from __future__ import annotations
import argparse,json
from pathlib import Path
from urllib.parse import quote,urljoin,urlparse
from urllib.request import urlopen,Request
API_ROOT="https://datadryad.org"
MAX_RESPONSE=2_000_000

def safe_dryad_api_url(url):
    p=urlparse(url)
    if p.scheme!="https" or p.netloc!="datadryad.org" or not p.path.startswith("/api/v2/"):
        raise ValueError("Only official Dryad v2 public metadata API supported")
    return url

def absolute_api_link(href):
    if not href or not isinstance(href,str):raise ValueError("Missing Dryad HAL link")
    return safe_dryad_api_url(urljoin(API_ROOT,href))

def api_json(url):
    req=Request(safe_dryad_api_url(url),headers={
       "Accept":"application/json",
       "User-Agent":"TTF-Quaternary-fresh-transplant-source-metadata-only-v01"})
    with urlopen(req,timeout=35) as resp:data=resp.read(MAX_RESPONSE+1)
    if len(data)>MAX_RESPONSE:raise ValueError("Dryad metadata exceeded size cap")
    obj=json.loads(data)
    if not isinstance(obj,dict):raise ValueError("Dryad metadata not an object")
    return obj

def select_version_link(dataset):
    return absolute_api_link(((dataset.get("_links") or {}).get("stash:version") or {}).get("href"))

def parse_one_dataset(dataset,doi,role):
    if dataset.get("identifier")!="doi:"+doi:
        raise ValueError("Dataset DOI differs from prospectively declared original")
    href=select_version_link(dataset)
    vid=urlparse(href).path.split("/")[-1]
    if not vid.isdigit():raise ValueError("Current Dryad version link not numeric")
    return {"doi":doi,"publication_role":role,"title":dataset.get("title"),
            "current_version_id":vid,"version_metadata_url":href,
            "publication_date":dataset.get("publicationDate"),
            "source_metadata_only":True,"phenotype_rows_inspected":False}

def parse_files_page(document):
    entries=(document.get("_embedded") or {}).get("stash:files")
    if not isinstance(entries,list):raise ValueError("Missing Dryad HAL stash:files")
    found=[]
    for e in entries:
        digest=e.get("digest")
        dtyp=str(e.get("digestType") or "").lower()
        dl=((e.get("_links") or {}).get("stash:download") or {}).get("href")
        if dl is not None:absolute_api_link(dl)
        found.append({
            "path":e.get("path"),"size_bytes":e.get("size"),
            "mime_type":e.get("mimeType"),"digest_type":e.get("digestType"),
            "catalog_digest":digest,
            "sha256_catalog":digest if dtyp in ("sha-256","sha256") else None,
            "download_href":dl,
            "data_records_read":False
        })
    nxt=(document.get("_links") or {}).get("next") or {}
    link=nxt.get("href") if isinstance(nxt,dict) else None
    return found,absolute_api_link(link) if link else None

def project_listing(candidate,dataset,files):
    result=parse_one_dataset(dataset,candidate["doi"],candidate["role"])
    match=[f for f in files
           if str(f.get("path") or "").lower().split("/")[-1]==candidate["expected_file_name"].lower()]
    result.update({
       "candidate_id":candidate["id"],"expected_file_name":candidate["expected_file_name"],
       "expected_file_present":bool(match),
       "source_catalog_status":"FILE_CATALOG_AND_EXPECTED_SOURCE_VERIFIED" if match else "EXPECTED_SOURCE_FILE_NOT_FOUND",
       "candidate_reason":candidate["why"],
       "source_overlap_with":candidate["known_overlap"],
       "independent_data_status":"UNKNOWN_POSSIBLE_OVERLAP",
       "file_count":len(files),"files":files,
       "raw_fitness_or_plasticity_values_read":False,
       "original_origin_coordinates_read":False,
       "cross_study_cohort_test_performed":False
    })
    return result

def read_catalog(candidate):
    encoded=quote("doi:"+candidate["doi"],safe="")
    d=api_json(f"{API_ROOT}/api/v2/datasets/{encoded}")
    ref=select_version_link(d)
    version=api_json(ref)
    filelink=((version.get("_links") or {}).get("stash:files") or {}).get("href")
    if not filelink:raise ValueError("Published version lacks file-list HAL link")
    url=absolute_api_link(filelink)
    found=[];seen=set()
    for _ in range(5):
        if url in seen:raise ValueError("Cyclic file catalog")
        seen.add(url)
        new,next_page=parse_files_page(api_json(url))
        found.extend(new)
        if not next_page:break
        url=next_page
    else:raise ValueError("More than five catalog pages: STOP")
    obj=project_listing(candidate,d,found)
    obj["version_api_file_listing_pages"]=len(seen)
    obj["version_source_status"]=version.get("versionStatus") or version.get("curationStatus")
    return obj

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--contract",required=True,type=Path)
    ap.add_argument("--output",required=True,type=Path)
    args=ap.parse_args()
    c=json.loads(args.contract.read_text())
    if c.get("status")!="FROZEN_BEFORE_EXTERNAL_DRYAD_METADATA_SEARCH":
        raise ValueError("No frozen fresh-source eligibility contract")
    results=[]
    for candidate in c["candidates"]:
        try:
            obj=read_catalog(candidate)
        except Exception as exc:
            obj={"candidate_id":candidate["id"],"doi":candidate["doi"],
                 "source_catalog_status":"SOURCE_METADATA_UNAVAILABLE_NOT_ECOLOGICAL_NULL",
                 "error_type":type(exc).__name__,"error_message":str(exc)[:220],
                 "independent_data_status":"UNKNOWN_POSSIBLE_OVERLAP",
                 "raw_fitness_or_plasticity_values_read":False,
                 "original_origin_coordinates_read":False}
        results.append(obj)
        print("DRYAD_SOURCE_METADATA_ONLY",candidate["id"],obj["source_catalog_status"],flush=True)
    summary={
        "schema":"ttf_fresh_reciprocal_transplant_dryad_catalog_v01",
        "status":"PUBLIC_SOURCE_CATALOG_CHECKED_NOT_BIOLOGICAL_EFFECT",
        "original_old_Haldane_terminal_still_closed":True,
        "source_coding_meta_only":True,
        "source_results":results,
        "raw_fitness_or_plasticity_values_read":False,
        "new_evolutionary_climate_predictive_effect_tested":False,
        "source_study_overlap_reconciled":False,
        "new_empirical_outcome_opening_authorized":False
    }
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(summary,sort_keys=True,ensure_ascii=False,indent=2)+"\n")
    print("CATALOG_COUNT",sum(z.get("expected_file_present",False) for z in results),"/",len(results))

if __name__=="__main__":main()
