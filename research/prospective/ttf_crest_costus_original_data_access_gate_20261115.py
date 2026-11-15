"""Costus 2022/2026 original Dryad source access and schema, no biology data opened.

Files are never committed or exported. For CSV, only the first header line is
decoded; for ZIP, only the archive filenames. SHA256/size must match the
previously audited public metadata. API 401/403 means NOT_EVALUABLE_ACCESS,
not zero phenotypes/pollinator observations.
"""
from __future__ import annotations
import csv
import hashlib
import io
import json
import urllib.error
import urllib.request
import zipfile

FILES={
 "costus_2022_observed_pollinators":{
    "file_id":1569887,"size":6994,
    "sha256":"e6bfa395a1792cfe648016609dbdaf8153170206d8ddc69f3e9f1023278562a5",
    "file":"Pollinator_observations.csv","content_kind":"csv_header",
 },
 "costus_2022_ancestry_tree":{
    "file_id":1569878,"size":90309,
    "sha256":"24f1cbfaa82ce4d33aad1ba80787cc611a6730aafb38f8244954c4d0d50ae7a1",
    "file":"chronogram_v8_2.tre","content_kind":"binary_checksum_only",
 },
 "costus_2026_efficiency_archive":{
    "file_id":4746470,"size":38955,
    "sha256":"00009bc6170473ddfc3e74003939a5c083c70e6d50268adfa9cc2a8279931796",
    "file":"Data.Availability(Dryad).zip","content_kind":"zip_member_names_only",
 },
}

def metadata_only(original:bytes,item:dict)->dict:
    expected=item["sha256"]
    actual=hashlib.sha256(original).hexdigest()
    if len(original)!=item["size"] or actual!=expected:
        raise ValueError("Original Dryad version file bytes drifted from SHA256/size")
    result={"status":"ORIGINAL_BYTES_VERIFIED_NO_OUTCOME_OPENED",
            "original_file":item["file"],"size":len(original),"sha256":actual}
    if item["content_kind"]=="csv_header":
        raw=original.decode("utf-8-sig")
        first=next(csv.reader(io.StringIO(raw)),None)
        if not first or len(first)!=len(set(first)):
            raise ValueError("missing or duplicated pollinator observation CSV header fields")
        result["header_only"]=first
        result["data_row_values_decoded"]=False
    elif item["content_kind"]=="zip_member_names_only":
        if not zipfile.is_zipfile(io.BytesIO(original)):
            raise ValueError("Original purported archive is not ZIP")
        with zipfile.ZipFile(io.BytesIO(original)) as z:
            result["archive_member_names"]=[x.filename for x in z.infolist()]
        result["archive_members_contents_decoded"]=False
    else:
        result["phylogenetic_tip_or_branch_values_read"]=False
    return result


def acquire(item:dict)->dict:
    url=f"https://datadryad.org/api/v2/files/{item['file_id']}/download"
    req=urllib.request.Request(url,headers={"User-Agent":"TTF-Costus-public-source-only/1.0"})
    try:
        with urllib.request.urlopen(req,timeout=35) as resp:
            raw=resp.read(item["size"]+1)
    except urllib.error.HTTPError as e:
        if e.code in (401,403,404):
            return {"status":"SOURCE_DOWNLOAD_NOT_ACCESSIBLE","http_status":e.code,
                    "official_api_url":url,"no_outcomes_opened":True,
                    "not_a_biological_non_detection":True}
        raise
    result=metadata_only(raw,item)
    result["official_api_url"]=url
    return result


def main():
    rows=[]
    for key,item in FILES.items():
        row=acquire(item)
        row.update({"source":key,"file_id":item["file_id"]})
        rows.append(row)
    print(json.dumps({"status":"COSTUS_SOURCE_ACCESS_ONLY",
                      "files":rows,"no_scientific_effects_computed":True,
                      "no_source_bytes_exported":True},indent=2,sort_keys=True))


if __name__=="__main__": main()
