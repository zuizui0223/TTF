"""Original PMC/bioRxiv Penstemon 2026 Supplement 1 PDF access and index.

The original Stone et al. DOI 10.64898/2026.07.21.739878 points to Table S1
in a supplemental PDF. This source-only routine checks a FINITE list of
candidate official publisher/archive URLs. It reads no voucher-table values,
no floral phenotypes, no SNP genotypes, no genomic effects.
"""
from __future__ import annotations
import hashlib
import io
import json
import urllib.error
import urllib.request
from pypdf import PdfReader

DOI="10.64898/2026.07.21.739878"
SOURCES={
    "PMC_supplement_generic":"https://pmc.ncbi.nlm.nih.gov/articles/PMC13419788/bin/media-1.pdf",
    "PMC_nihpp_supplement":"https://pmc.ncbi.nlm.nih.gov/articles/PMC13419788/bin/nihpp-2026.07.21.739878v1-supp1.pdf",
    "bioRxiv_original_supplement":"https://www.biorxiv.org/content/10.64898/2026.07.21.739878v1.supplementary-material",
    "bioRxiv_early_supplement":"https://www.biorxiv.org/content/biorxiv/early/2026/07/24/2026.07.21.739878/DC1/embed/media-1.pdf?download=1",
}
MAX_PDF_BYTES=15_000_000
MAX_PDF_PAGES=80


def inspect_pdf_bytes(raw:bytes)->dict:
    if not raw.startswith(b"%PDF-"):
        raise ValueError("official source response is not a PDF signature")
    if len(raw)>MAX_PDF_BYTES:raise ValueError("original source PDF exceeds length cap")
    r=PdfReader(io.BytesIO(raw),strict=False)
    n=len(r.pages)
    if n<1 or n>MAX_PDF_PAGES:
        raise ValueError("original supplemental PDF page count outside source safety limit")
    # Qualify the table index only. Never transcribe species/voucher row values.
    headings={"Table S1":[],"Table S2":[],"Table S3":[]}
    for i,page in enumerate(r.pages):
        text=page.extract_text() or ""
        for heading in headings:
            if heading.casefold() in text.casefold():
                headings[heading].append(i+1)
    return {
       "status":"ORIGINAL_SUPPLEMENT_PDF_CONTENT_VERIFIED_INDEX_ONLY",
       "sha256":hashlib.sha256(raw).hexdigest(),
       "bytes":len(raw),
       "page_count":n,
       "table_index_mentions_page_number_1_based":headings,
       "no_original_specimen_rows_or_phenotype_values_exported":True,
       "source_pdf_not_equivalent_to_34_phenotype_specimen_manifest":True,
    }


def fetch_candidate(url):
    request=urllib.request.Request(url,headers={
       "Accept":"application/pdf",
       "User-Agent":"TTF-Penstemon-source-table-metadata/1.0",
    })
    try:
        with urllib.request.urlopen(request,timeout=40) as resp:
            data=resp.read(MAX_PDF_BYTES+1)
            content_type=resp.headers.get("Content-Type","")
    except urllib.error.HTTPError as e:
        return {"status":"OFFICIAL_SOURCE_HTTP_BLOCKED","http_status":e.code}
    except urllib.error.URLError as e:
        return {"status":"OFFICIAL_SOURCE_CONNECTION_UNAVAILABLE",
                "error_type":type(e).__name__}
    if len(data)>MAX_PDF_BYTES:
        return {"status":"SOURCE_BYTES_UNEXPECTEDLY_LARGE"}
    if not data.startswith(b"%PDF-"):
        return {"status":"OFFICIAL_SOURCE_RESPONSE_NOT_PDF",
                "content_type":content_type,"bytes":len(data)}
    return inspect_pdf_bytes(data)


def main():
    reports={name:fetch_candidate(url) for name,url in SOURCES.items()}
    verified=[k for k,v in reports.items() if v["status"]=="ORIGINAL_SUPPLEMENT_PDF_CONTENT_VERIFIED_INDEX_ONLY"]
    print(json.dumps({
       "status":"PENSTEMON_ORIGINAL_SUPPLEMENT_SOURCE_METADATA_ONLY",
       "paper_doi":DOI,
       "candidate_source_statuses":reports,
       "accessible_official_supplement_candidates":verified,
       "phenotyped_34_source_specimen_list_extracted":False,
       "no_floral_or_genomic_outcome_opened":True,
    },indent=2,sort_keys=True))


if __name__=="__main__":
    main()
