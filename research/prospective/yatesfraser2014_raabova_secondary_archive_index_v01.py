"""DOI-exact secondary published Raabova-data archive discovery, metadata only.

Yates & Fraser (2014) explicitly acknowledge J. Raabova as a provider of data
and archive their meta-analysis to Dryad 10.5061/dryad.81k4v. This is NOT the
original 2007 adult transplant plant-tag table. Query only Dryad metadata and
the official Zenodo Dryad preservation index, never read archive file bodies.
"""
from __future__ import annotations
import json,re
from pathlib import Path
from urllib.error import HTTPError,URLError
from urllib.parse import urlencode
from urllib.request import Request,urlopen

CONTRACT=Path(__file__).with_name("yatesfraser2014_raabova_secondary_dryad_index_contract_v01.json")

def doi(value:object)->str:
    if not isinstance(value,str):return ""
    text=value.strip().lower()
    return re.sub(r"^(?:https?://(?:dx\.)?doi\.org/|doi:)", "",text).rstrip("/")

def exact_original_doi(record:dict,target:str)->bool:
    meta=record.get("metadata") or {}
    v=[]
    for key in ("related_identifiers","identifiers","alternate_identifiers"):
        for identifier in meta.get(key) or []:
            if isinstance(identifier,dict):
                v.append(identifier.get("identifier"))
    v.extend([record.get("doi"),meta.get("doi")])
    pids=record.get("pids") or {}
    if isinstance(pids,dict) and isinstance(pids.get("doi"),dict):
        v.append(pids["doi"].get("identifier"))
    return target in {doi(x) for x in v}

def safe_files(record:dict)->list[dict]:
    raw=record.get("files") or []
    if isinstance(raw,dict):
        items=list(raw.get("entries",{}).values()) if isinstance(raw.get("entries"),dict) else []
    elif isinstance(raw,list):
        items=raw
    else:
        items=[]
    files=[]
    for item in items:
        if not isinstance(item,dict):continue
        name=item.get("key") or item.get("filename")
        if not isinstance(name,str) or len(name)>260:continue
        files.append({"filename":name,"size_bytes":item.get("size"),
             "published_checksum":item.get("checksum")})
    return sorted(files,key=lambda x:x["filename"])

def zenodo_meta(c:dict)->dict:
    conf=c["sources"];base=conf["zenodo_official_records_api"]
    target=c["origin"]["archived_dataset_doi"]
    found={};attempts=[]
    for phrase in conf["zenodo_exact_doi_queries"]:
        endpoint=base+"?"+urlencode({"q":phrase,"size":conf["max_zenodo_search_records"]})
        try:
            with urlopen(Request(endpoint,headers={
                "Accept":"application/json",
                "User-Agent":"TTF-Raabova-original-ID-secondary-Dryad-preservation-metadata-only/2026-11-16"
            }),timeout=30) as resp:
                raw=resp.read(conf["max_source_metadata_bytes"]+1)
            if len(raw)>conf["max_source_metadata_bytes"]:raise ValueError("Zenodo metadata response cap exceeded")
            hits=json.loads(raw).get("hits",{}).get("hits",[])
            if not isinstance(hits,list):raise ValueError("Zenodo metadata hits schema")
            attempts.append({"query":phrase,"status":"SEARCH_RETURNED","hits_reviewed":len(hits)})
            for hit in hits:
                if not isinstance(hit,dict) or not exact_original_doi(hit,target):continue
                rid=hit.get("id")
                if not isinstance(rid,int) or rid<1:continue
                found[rid]={"id":rid,"url":f"https://zenodo.org/records/{rid}",
                     "original_Dryad_DOI_exact_verified":True,
                     "metadata_file_list":safe_files(hit)}
        except (HTTPError,URLError,ValueError,TimeoutError,json.JSONDecodeError) as exc:
            attempts.append({"query":phrase,"status":"ACCESS_OR_METADATA_ERROR",
                 "error_type":type(exc).__name__,"http_status":getattr(exc,"code",None)})
    return {"matched_exact_original_Dryad_doi_preservation_records":
            [found[k] for k in sorted(found)],"attempts":attempts}

def dryad_meta(c:dict)->dict:
    endpoint=c["sources"]["dryad_official_api"]
    maxbytes=c["sources"]["max_source_metadata_bytes"]
    try:
        with urlopen(Request(endpoint,headers={
            "Accept":"application/json","User-Agent":"TTF-publisher-Dryad-original-file-list-only/2026-11-16"
        }),timeout=30) as response:
            raw=response.read(maxbytes+1)
        if len(raw)>maxbytes:raise ValueError("Dryad metadata cap")
        x=json.loads(raw)
        found_doi=x.get("identifier") or x.get("doi")
        if doi(found_doi)!=c["origin"]["archived_dataset_doi"]:
            raise ValueError("Dryad returned different original dataset DOI")
        return {
          "status":"EXACT_DATASET_METADATA_FOUND",
          "source_doi_exact_verified":True,
          "dataset_title":str(x.get("title") or "")[:200],
          "original_file_links_NOT_downloaded":True,
          "public_metadata_file_list":safe_files(x)
        }
    except (HTTPError,URLError,TimeoutError,ValueError,json.JSONDecodeError) as exc:
        return {"status":"OFFICIAL_DRYAD_API_METADATA_ACCESS_NOT_EVALUABLE",
            "error_type":type(exc).__name__,"http_status":getattr(exc,"code",None),
            "record_absence_not_inferred":True}

def probe(c:dict)->dict:
    if c.get("status")!="FROZEN_BEFORE_SECONDARY_ARCHIVE_SOURCE_FILE_METADATA_SEARCH":
        raise ValueError("No precommitted archive-index access rule")
    if c["sources"]["data_file_download_authorized"] or not c["exact_Dryad_DOI_required_before_accepted_zenodo_candidate"]:
        raise ValueError("Data files and approximate DOI matching are not authorized")
    z=zenodo_meta(c)
    d=dryad_meta(c)
    return {
        "schema":"ttf_YatesFraser2014_Raabova_shared_secondary_dataset_DOI_exact_metadata_preflight_v01",
        "source_paper_doi":c["origin"]["paper_doi"],
        "archived_secondary_Dryad_doi":c["origin"]["archived_dataset_doi"],
        "paper_acknowledges_Raabova_input_data":True,
        "origin_published_analysis_874_estimates_across_13_species_not_2007_individual_ID_rows":True,
        "zenodo_preservation":z,"dryad_original_metadata":d,
        "status":("EXACT_DOI_SOURCE_INDEX_PRESENT_METADATA_ONLY"
            if z["matched_exact_original_Dryad_doi_preservation_records"] or
               d["status"]=="EXACT_DATASET_METADATA_FOUND"
            else "SOURCE_METADATA_NO_CONFIRMED_INDEX_NOT_PROOF_OF_ABSENCE"),
        "raw_data_files_or_spreadsheets_downloaded":False,
        "Raabova2007_original_field_individual_mark_column_verified":False,
        "2007_vs_2014_individual_overlap_count":None,
        "independent_new_experiment_certified":False,
        "old_bontrager_44of50_and_2014_source_only_lineage_unchanged":True,
        "fitness_genetic_LGM_adaptation_plasticity_speciation_outcomes_opened":False
    }

if __name__=="__main__":
    c=json.loads(CONTRACT.read_text(encoding="utf-8"))
    print("YATESFRASER_RAABOVA_2014_SECONDARY_ARCHIVE_INDEX="+json.dumps(probe(c),sort_keys=True))
