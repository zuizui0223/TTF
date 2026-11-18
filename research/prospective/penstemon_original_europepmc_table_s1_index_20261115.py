"""Authentic Europe PMC supplemental ZIP: Penstemon Table S1 PDF header index.

Original preprint PMC13419788, public official Europe PMC supplementaryFiles.
The archive SHA256 was verified in source-only run 43711613364. Extract only
PDF byte identity and Table S1/S2/S3 page headings, not sample rows,
not floral measurements, not genetic outcomes.
"""
from __future__ import annotations
import hashlib
import io
import json
import re
import urllib.request
import zipfile
from pypdf import PdfReader

URL="https://www.ebi.ac.uk/europepmc/webservices/rest/PMC13419788/supplementaryFiles"
# Europe PMC dynamically rebuilds ZIP wrappers: prior archive digest is
# provenance for an earlier manifest, NOT a stable PDF member digest.
EARLIER_ZIP_SHA="42abf929e1428e9f643e76ed80057574b66c4bcd2c47b77ff50ef370c4a01278"
PDF_MEMBER="media-1.pdf"
PDF_DECLARED_BYTES=2808017
PDF_MAX_BYTES=7_000_000
SOURCE_ZIP_MAX=6_000_000


def index_pdf(raw:bytes)->dict:
    if not raw.startswith(b"%PDF-") or len(raw)!=PDF_DECLARED_BYTES:
        raise ValueError("not original publisher supplemental PDF bytes")
    reader=PdfReader(io.BytesIO(raw),strict=False)
    pages={}
    for i,page in enumerate(reader.pages,1):
        tx=page.extract_text() or ""
        events=[]
        for label in ["Table S1","Table S2","Table S3","Table S4","Supplementary Table 1","Phenotyping","Sampling"]:
            if re.search(r"\b"+re.escape(label)+r"\b",tx,re.I):
                events.append(label)
        if events:
            # Only labeled page hits, not numerical outcome or specimen rows.
            pages[str(i)]=events
    return {"status":"OFFICIAL_PDF_PAGES_INDEXED_NO_TABLE_ROWS",
            "pdf_sha256":hashlib.sha256(raw).hexdigest(),
            "pdf_bytes":len(raw),"page_count":len(reader.pages),
            "pages_with_author_table_headings":pages,
            "no_specimen_or_flower_value_export":True}


EXPECTED_MEMBERS={
    "media-1.pdf",
    "nihpp-2026.07.21.739878v1-f0001.gif",
    "nihpp-2026.07.21.739878v1-f0002.gif",
    "nihpp-2026.07.21.739878v1-f0003.gif",
    "nihpp-2026.07.21.739878v1-f0004.gif",
    "nihpp-2026.07.21.739878v1-f0001.jpg",
    "nihpp-2026.07.21.739878v1-f0002.jpg",
    "nihpp-2026.07.21.739878v1-f0003.jpg",
    "nihpp-2026.07.21.739878v1-f0004.jpg",
}


def extract_original_pdf_member(archive:bytes) -> bytes:
    if len(archive)>SOURCE_ZIP_MAX or not zipfile.is_zipfile(io.BytesIO(archive)):
        raise ValueError("Europe PMC ZIP unavailable or unsafe")
    with zipfile.ZipFile(io.BytesIO(archive)) as z:
        members={info.filename:info for info in z.infolist() if not info.is_dir()}
        if set(members)!=EXPECTED_MEMBERS:
            raise ValueError("Europe PMC original supplemental ZIP member inventory drift")
        if members[PDF_MEMBER].file_size!=PDF_DECLARED_BYTES:
            raise ValueError("original PDF member length differs from earlier official manifest")
        with z.open(PDF_MEMBER) as source:
            raw=source.read(PDF_MAX_BYTES+1)
    if len(raw)!=PDF_DECLARED_BYTES:
        raise ValueError("original PDF uncompressed byte length drift")
    if not raw.startswith(b"%PDF-"):
        raise ValueError("original PDF member is not PDF")
    return raw


def source_get():
    req=urllib.request.Request(URL,headers={
       "User-Agent":"TTF-Habroanthus-EPMC-original-PDF-index/1.0",
       "Cache-Control":"no-cache",
    })
    with urllib.request.urlopen(req,timeout=45) as response:
        data=response.read(SOURCE_ZIP_MAX+1)
    return data


def main():
    a=source_get()
    b=source_get()
    pdf_a=extract_original_pdf_member(a)
    pdf_b=extract_original_pdf_member(b)
    hash_a=hashlib.sha256(pdf_a).hexdigest()
    hash_b=hashlib.sha256(pdf_b).hexdigest()
    if hash_a!=hash_b:
        raise ValueError("official PDF member hash differs between two independent requests")
    result=index_pdf(pdf_a)
    print(json.dumps({
       "status":"PENSTEMON_OFFICIAL_SUPPLEMENT_TABLE_INDEX_ONLY",
       "original_source_url":URL,
       "archive_wrapper_sha256_per_request":[hashlib.sha256(a).hexdigest(),
                                             hashlib.sha256(b).hexdigest()],
       "earlier_source_zip_sha256":EARLIER_ZIP_SHA,
       "archive_wrapper_byte_equality":a==b,
       "pdf_member_same_sha256_in_two_requests":hash_a==hash_b,
       "source_member":PDF_MEMBER,
       "pdf_metadata":result,
       "paper_reported_sequenced_accessions":52,
       "paper_reported_flower_phenotyped_accessions":34,
       "identified_34_phenotyped_BioSample_IDs":False,
       "no_floral_or_genetic_outcomes_read":True,
    },indent=2,sort_keys=True))


if __name__=="__main__":
    main()
