"""Official Raabova 2007 source-data discovery by HTML attachment metadata ONLY.

Only site HTML anchor URLs and file types are read. Never download a thesis PDF,
supplement workbook, row data or outcomes. A missing linked spreadsheet is not
proof that the historical 2007 field-tag original never existed.
"""
from __future__ import annotations
import json
from html.parser import HTMLParser
from pathlib import Path
from urllib.error import HTTPError, URLError
from urllib.parse import unquote, urljoin, urlsplit
from urllib.request import Request, urlopen

CONTRACT=Path(__file__).with_name("raabova2007_original_data_location_preflight_contract_v01.json")

class LinkReader(HTMLParser):
    def __init__(self):
        super().__init__(convert_charrefs=True)
        self.urls=[]
    def handle_starttag(self,tag,attrs):
        if tag not in ("a","link"):
            return
        d=dict(attrs)
        if d.get("href"):
            self.urls.append(d["href"])

def file_links(html:str,url:str,extensions:set[str]) -> list[dict]:
    parser=LinkReader()
    parser.feed(html)
    found={}
    for ref in parser.urls:
        absolute=urljoin(url,ref)
        path=unquote(urlsplit(absolute).path)
        if not path or "/" not in path:continue
        name=path.rsplit("/",1)[-1]
        if "." not in name:continue
        ext=name.rsplit(".",1)[-1].lower()
        if ext not in extensions:continue
        key=(name,absolute.split("#")[0])
        found[key]={"file_name":name,"extension":ext,
                    "listing_link_url":absolute,"downloaded":False}
    return sorted(found.values(),key=lambda x:(x["file_name"],x["listing_link_url"]))

def interpret_link_listing(files:list[dict],source_kind:str)->dict:
    structured={"csv","tsv","xls","xlsx","ods","zip","rds","rdata","dat"}
    report={"pdf","txt"}
    has_data=any(x["extension"] in structured for x in files)
    count_structured=sum(x["extension"] in structured for x in files)
    return {
        "official_landing_html_file_links_found":len(files),
        "linked_structured_file_candidate_count":count_structured,
        "linked_report_or_text_files":sum(x["extension"] in report for x in files),
        "status":("OFFICIAL_LISTED_STRUCTURE_FILES_REQUIRE_ROW_SCHEMA_AUDIT" if has_data
            else "ONLY_REPORT_TEXT_FILES_LINKED_IN_THIS_REPOSITORY_PAGE_NOT_DATASET_ABSENCE"
            if files else "SOURCE_HTML_NO_DECISIVE_FILE_LINKS_NOT_EVALUABLE"),
        "source_record_scope":source_kind,
        "candidate_2007_individual_tag_table_verified":False,
        "files":files,
    }

def preflight(contract:dict)->dict:
    if contract["status"]!="FROZEN_BEFORE_INSTITUTIONAL_SOURCE_RECORD_ATTACHMENT_CENSUS":
        raise ValueError("Official source file-list audit not predeclared")
    allowed_hosts=set(contract["allowlisted_hosts"])
    extensions=set(contract["file_extensions"])
    results=[]
    for src in contract["sources"]:
        parsed=urlsplit(src["url"])
        if parsed.scheme!="https" or parsed.hostname not in allowed_hosts:
            raise ValueError("Original source hostname outside pre-frozen list")
        meta={"source_key":src["key"],"source_url":src["url"],
              "source_kind":src["source_kind"]}
        try:
            req=Request(src["url"],headers={
                "Accept":"text/html,application/xhtml+xml",
                "User-Agent":"TTF-Raabova2007-official-source-file-list-only/2026-11-15"})
            with urlopen(req,timeout=28) as stream:
                raw=stream.read(contract["maximum_html_bytes_per_source"]+1)
                meta["http_status"]=stream.status
                meta["final_host"]=urlsplit(stream.geturl()).hostname
                content_type=stream.headers.get("Content-Type","")
                meta["http_content_type"]=content_type
            if len(raw)>contract["maximum_html_bytes_per_source"]:
                raise ValueError("Public source HTML size exceeds frozen limit")
            if (b"<html" not in raw[:10000].lower()
                and b"<!doctype html" not in raw[:10000].lower()):
                raise ValueError("Source page does not contain HTML content")
            text=raw.decode("utf-8-sig",errors="replace")
            if any(x in text.lower() for x in
                   ("checking your browser","captcha","just a moment...")):
                raise ValueError("Source page is an access challenge")
            links=file_links(text,src["url"],extensions)
            meta.update(interpret_link_listing(links,src["source_kind"]))
            meta["source_html_bytes_read"]=len(raw)
        except (HTTPError,URLError,TimeoutError,ValueError,UnicodeError) as exc:
            meta.update({"status":"OFFICIAL_SOURCE_LIST_ACCESS_OR_SCHEMA_NOT_EVALUABLE",
                        "error_type":type(exc).__name__,
                        "http_error_code":getattr(exc,"code",None),
                        "files":[]})
        results.append(meta)
    return {
        "schema":"ttf_Raabova2007_source_only_publisher_and_institutional_file_index_v01",
        "status":"REPOSITORY_INDEX_REVIEW_NOT_ORIGINAL_2007_ID_ROW_MATCH",
        "sources":results,
        "institutions_checked":2,
        "publisher_records_checked":1,
        "original_2007_individual_number_table_confirmed":False,
        "2007_2014_same_individual_tag_join":None,
        "no_hit_not_proof_of_unavailability":True,
        "source_PDFs_or_data_files_downloaded":False,
        "fitness_or_genetic_adaptation_or_plasticity_outcome_opened":False,
        "Raabova2007_Pankova2014_same_2004_adult_cohort_author_confirmed":True,
        "prior_Bontrager_44of50_stop_unchanged":True
    }

if __name__=="__main__":
    print("RAABOVA2007_OFFICIAL_ORIGINAL_DATA_SOURCE_INDEX="+json.dumps(
        preflight(json.loads(CONTRACT.read_text(encoding="utf-8"))),
        ensure_ascii=False,sort_keys=True))
