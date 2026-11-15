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
ZIP_SHA="42abf929e1428e9f643e76ed80057574b66c4bcd2c47b77ff50ef370c4a01278"
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


def main():
    req=urllib.request.Request(URL,headers={"User-Agent":"TTF-Habroanthus-EPMC-original-PDF-index/1.0"})
    with urllib.request.urlopen(req,timeout=45) as response:
        archive=response.read(SOURCE_ZIP_MAX+1)
    if len(archive)>SOURCE_ZIP_MAX or hashlib.sha256(archive).hexdigest()!=ZIP_SHA:
        raise ValueError("Europe PMC original supplementary ZIP SHA256 mismatch")
    with zipfile.ZipFile(io.BytesIO(archive)) as z:
        members={info.filename:info for info in z.infolist() if not info.is_dir()}
        if PDF_MEMBER not in members:
            raise ValueError("expected original supplemental PDF member missing")
        if members[PDF_MEMBER].file_size!=PDF_DECLARED_BYTES:
            raise ValueError("original PDF member declared length drift")
        with z.open(PDF_MEMBER) as member:
            pdf_bytes=member.read(PDF_MAX_BYTES+1)
    if len(pdf_bytes)>PDF_MAX_BYTES:raise ValueError("original PDF over cap")
    result=index_pdf(pdf_bytes)
    print(json.dumps({
       "status":"PENSTEMON_OFFICIAL_SUPPLEMENT_TABLE_INDEX_ONLY",
       "original_source_url":URL,
       "archive_sha256":ZIP_SHA,
       "source_member":PDF_MEMBER,
       "pdf_metadata":result,
       "paper_reported_sequenced_accessions":52,
       "paper_reported_flower_phenotyped_accessions":34,
       "identified_34_phenotyped_BioSample_IDs":False,
       "no_floral_or_genetic_outcomes_read":True,
    },indent=2,sort_keys=True))


if __name__=="__main__":main()
