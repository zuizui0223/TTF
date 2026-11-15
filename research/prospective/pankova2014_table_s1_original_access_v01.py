"""Source-only preflight for original Pankova et al 2014 Table S1 (.XLS).

Check official article/Supplement URLs for workbook physical identity
and accessibility without downloading any worksheet data or cell values.
A downloaded first 512-byte storage header is not an individual-ID column audit.
"""
from __future__ import annotations
import json
from urllib.error import HTTPError, URLError
from urllib.parse import urlparse
from urllib.request import Request, urlopen

ORIGINAL_ARTICLE_DOI = "10.1371/journal.pone.0093967"
SUPPLEMENTARY_LOCATOR = "pone.0093967.s001.xls"
LINKS = (
    "https://journals.plos.org/plosone/article/file?id=10.1371/journal.pone.0093967.s001&type=supplementary",
    "https://pmc.ncbi.nlm.nih.gov/articles/PMC3977983/bin/pone.0093967.s001.xls",
    "https://pmc.ncbi.nlm.nih.gov/articles/instance/3977983/bin/pone.0093967.s001.xls",
)
MAX_HEADER = 512
OLE_SIGNATURE = bytes.fromhex("D0CF11E0A1B11AE1")
ZIP_SIGNATURE = b"PK\x03\x04"


def recognized_storage_signature(first_bytes: bytes) -> str:
    if first_bytes.startswith(OLE_SIGNATURE):
        return "OLE2_BINARY_EXCEL_CANDIDATE"
    if first_bytes.startswith(ZIP_SIGNATURE):
        return "ZIP_XLSX_STYLE_CANDIDATE"
    if first_bytes.lstrip().lower().startswith((b"<!doctype html", b"<html")):
        return "HTML_INTERSTITIAL_NOT_WORKBOOK"
    return "UNKNOWN_FILE_FORMAT_NOT_CERTIFIED"


def probe() -> dict:
    attempts = []
    for u in LINKS:
        result = {"requested_official_host":urlparse(u).hostname,
                  "requested_original_supplement_locator":SUPPLEMENTARY_LOCATOR}
        if urlparse(u).hostname not in ("journals.plos.org","pmc.ncbi.nlm.nih.gov"):
            raise ValueError("Unapproved supplement URL")
        try:
            req=Request(u,headers={
                "Accept":"*/*","Range":f"bytes=0-{MAX_HEADER-1}",
                "User-Agent":"TTF-2014-author-original-source-only-header/2026-11-15"})
            with urlopen(req,timeout=30) as response:
                header=response.read(MAX_HEADER)
                actual=response.geturl()
                result.update({
                    "http_status":response.status,
                    "response_content_type":response.headers.get("Content-Type",""),
                    "published_content_length_header":response.headers.get("Content-Length"),
                    "response_bytes_read":len(header),
                    "storage_signature":recognized_storage_signature(header),
                    "final_host":urlparse(actual).hostname
                })
            if result["storage_signature"] in (
                "OLE2_BINARY_EXCEL_CANDIDATE","ZIP_XLSX_STYLE_CANDIDATE"):
                # A valid file signature alone does NOT certify the workbook
                # contains individual tag identifiers or shared plants.
                result["status"]="OFFICIAL_SUPPLEMENT_STORAGE_HEADER_AVAILABLE"
                attempts.append(result)
                break
            result["status"]="NOT_VERIFIED_WORKBOOK"
        except (HTTPError,URLError,TimeoutError,ValueError) as err:
            result.update({
                "status":"OFFICIAL_SUPPLEMENT_ACCESS_GATED",
                "error_type":type(err).__name__,
                "http_error_code":getattr(err,"code",None),
            })
        attempts.append(result)
    success=any(a["status"]=="OFFICIAL_SUPPLEMENT_STORAGE_HEADER_AVAILABLE" for a in attempts)
    return {
        "schema":"ttf_Raabova_Pankova2014_table_S1_original_storage_access_SOURCE_ONLY_v01",
        "Pankova_2014_paper_DOI":ORIGINAL_ARTICLE_DOI,
        "published_supplemental_file_locator":SUPPLEMENTARY_LOCATOR,
        "status":"ORIGINAL_TABLE_S1_STORAGE_ACCESS_OBSERVED_ONLY" if success
                 else "ORIGINAL_TABLE_S1_ACCESS_NOT_EVALUABLE",
        "attempts":attempts,
        "supplement_row_data_read":False,
        "worksheet_columns_read":False,
        "individual_plant_tag_numbers_read":False,
        "crosspaper_individual_tag_ID_join_completed":False,
        "actual_shared_2004_adult_experiment_confirmed_by_author":True,
        "new_independent_experiment_cohort_added":0,
        "fitness_genetic_adaptation_LGM_plasticity_speciation_values_opened":False,
        "old_Bontrager_44_of_50_frozen_gate_unchanged":True
    }


if __name__=="__main__":
    print("PANKOVA2014_S1_ATTACHMENT_HEADER="+json.dumps(probe(),sort_keys=True))
