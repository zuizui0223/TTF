"""Costus independent pollinator-transition source discovery: metadata only.

2022 ancestral-tree/pollinator/morphology Dryad source, 2025 scent OSF node,
and 2026 independent pollinator-efficiency Dryad source. This intentionally
does NOT download phenotypes, visitation responses, genetic sequences or scents.
"""
from __future__ import annotations
import json
import urllib.error
import urllib.parse
import urllib.request

DATASETS = {
 "kay_grossenbacher_2022":{
    "doi":"10.7291/D1C39G",
    "expected_files":{
        "chronogram_v8_2.tre","Pollinator_observations.csv",
        "Costus_traits_greenhouse.csv","Costus_traits_Maas.csv",
    },
 },
 "costus_efficiency_2026":{
    "doi":"10.5061/dryad.931zcrk0r",
    "expected_files":set(), # Archive names are not inferred before source inspection.
 },
}
OSF_NODE="2ap4k"
BASE="https://datadryad.org/api/v2"


def fetch_json(url: str):
    req=urllib.request.Request(url,headers={
      "Accept":"application/json","User-Agent":"TTF-CREST-Costus-metadata-audit/1.0",
    })
    try:
        with urllib.request.urlopen(req, timeout=45) as response:
            data=response.read(2_000_001)
    except urllib.error.HTTPError as err:
        return {"status":"NOT_ACCESSIBLE","http_status":err.code,"url":url}
    if len(data)>2_000_000:raise ValueError("metadata JSON over length cap")
    return json.loads(data)


def href(x):
    if isinstance(x,str):return x
    if isinstance(x,dict):return x.get("href","")
    return ""


def hal_items(payload:dict,term:str) -> list:
    embedded=payload.get("_embedded",{})
    for key,val in embedded.items():
        if term in key.lower() and isinstance(val,list):
            return val
    raise ValueError(f"HAL collection {term} not found in {sorted(embedded)}")


def parse_versions(payload:dict)->list[dict]:
    rows=[]
    for v in hal_items(payload,"versions"):
        link=href(v.get("_links",{}).get("self"))
        ending=link.rsplit("/",1)[-1]
        version_id=v.get("id") or (int(ending) if ending.isdigit() else None)
        if not isinstance(version_id,int):raise ValueError("Dryad version ID unresolved")
        rows.append({
          "id":version_id,"version_number":v.get("versionNumber"),
          "publication_date":v.get("publicationDate"),
          "self":link,
        })
    if not rows:raise ValueError("zero source versions")
    return rows


def parse_files(payload:dict)->list[dict]:
    rows=[]
    for f in hal_items(payload,"files"):
        link=href(f.get("_links",{}).get("self"))
        ending=link.rsplit("/",1)[-1]
        file_id=f.get("id") or (int(ending) if ending.isdigit() else None)
        rows.append({
           "file_id":file_id,"filename":f.get("fileName") or f.get("path"),
           "bytes":f.get("size"),"digest_type":f.get("digestType"),
           "digest":f.get("digest"),
           "self":link,
        })
    return rows


def dryad_metadata(key:str,item:dict)->dict:
    doi=item["doi"]
    encoded=urllib.parse.quote("doi:"+doi,safe="")
    versions_url=f"{BASE}/datasets/{encoded}/versions"
    a=fetch_json(versions_url)
    if a.get("status")=="NOT_ACCESSIBLE":
        return {"source":key, "doi":doi, "status":"VERSION_METADATA_NOT_ACCESSIBLE",
                "http_status":a.get("http_status")}
    versions=parse_versions(a)
    sel=max(versions,key=lambda x:(int(x["version_number"] or 0),x["id"]))
    b=fetch_json(f"{BASE}/versions/{sel['id']}/files")
    if b.get("status")=="NOT_ACCESSIBLE":
        return {"source":key,"doi":doi,"status":"FILE_METADATA_NOT_ACCESSIBLE",
                "selected_version":sel,"http_status":b.get("http_status")}
    files=parse_files(b)
    found={f["filename"] for f in files}
    missing=sorted(item["expected_files"]-found)
    return {
       "source":key,"doi":doi,"status":"SOURCE_FILE_METADATA_CHECKED" if not missing
         else "EXPECTED_SOURCE_FILES_MISSING_HOLD",
       "versions":versions, "selected_version":sel, "files":files,
       "expected_source_file_missing":missing,"raw_files_downloaded":False,
    }


def osf_metadata()->dict:
    url=f"https://api.osf.io/v2/nodes/{OSF_NODE}/"
    response=fetch_json(url)
    if response.get("status")=="NOT_ACCESSIBLE":
        return {"source":"darragh_2025_scent_osf",
                "status":"OSF_METADATA_NOT_ACCESSIBLE","http_status":response["http_status"],
                "node_id":OSF_NODE,"node_url":f"https://osf.io/{OSF_NODE}/"}
    data=response.get("data",{})
    if data.get("id")!=OSF_NODE:
        raise ValueError("OSF returned unexpected project id")
    attrs=data.get("attributes",{})
    return {"source":"darragh_2025_scent_osf","status":"OSF_PROJECT_METADATA_ONLY",
            "node_id":data["id"],"title":attrs.get("title"),
            "public":attrs.get("public"),"node_url":f"https://osf.io/{OSF_NODE}/",
            "raw_files_downloaded":False}


def main():
    out=[dryad_metadata(k,v) for k,v in DATASETS.items()]
    out.append(osf_metadata())
    print(json.dumps({
      "status":"COSTUS_PRIMARY_SOURCE_METADATA_ONLY_NO_OUTCOMES",
      "sources":out,
      "no_species_crosswalk_performed":True,
      "no_observed_pollinator_outcomes_read":True,
      "no_scent_or_molecular_responses_read":True,
      "cross_study_independence_unverified":True,
    },indent=2,sort_keys=True,ensure_ascii=False))


if __name__=="__main__":main()
