"""Original 2026 Penstemon Table S1 source specimen-header diagnostic.

Reads only page 11 of SHA256-pinned media-1.pdf distributed inside the official
Europe PMC supplemental ZIP. Does not read floral phenotype trait values,
genotypes, evolutionary effects, or assume source voucher-to-SAMN mapping.
"""
from __future__ import annotations
import hashlib
import io
import json
import re
import sys
from pathlib import Path
from pypdf import PdfReader

HERE=str(Path(__file__).resolve().parent)
if HERE not in sys.path:sys.path.insert(0,HERE)
from penstemon_original_europepmc_table_s1_index_20261115 import (
    source_get,extract_original_pdf_member
)

ORIGINAL_PDF_SHA="3682720dad4f0621bebef3be7014aa667067c5f0dce9160ca8f4c4a758e7f943"
TABLE_S1_PAGE_1BASED=11


def inspect_s1_layout(text:str)->dict:
    if not re.search(r"\bTable\s*S1\b",text,re.I):
        raise ValueError("original Table S1 page heading missing")
    nonempty=[x.strip() for x in text.splitlines() if x.strip()]
    if not nonempty:
        raise ValueError("no selectable Table S1 source text")
    # Short bounded original source schema preview only. The full 52-row table
    # is never exported to the repo or logs.
    preview=nonempty[:20]
    return {
        "status":"ORIGINAL_TABLE_S1_HEADER_AND_LAYOUT_PREVIEW_ONLY",
        "source_nonempty_lines":len(nonempty),
        "source_first_20_lines":preview,
        "sample_author_key_prefix_occurrences":{
           "CAW":len(re.findall(r"\bCAW[- ]\d{2,4}\b",text)),
           "BWS":len(re.findall(r"\bBWS[- ]\d{2,4}\b",text)),
           "LW":len(re.findall(r"\bLW[- ]\d{2,4}\b",text)),
           "SAMN":len(re.findall(r"\bSAMN\d+\b",text)),
        },
        "source_table_proves_34_phenotyped":"NOT_ASSESSED",
        "no_floral_trait_or_genotype_values_read":True,
    }


def main():
    raw=extract_original_pdf_member(source_get())
    if hashlib.sha256(raw).hexdigest()!=ORIGINAL_PDF_SHA:
        raise ValueError("original 2026 supplemental media-1.pdf SHA256 differs")
    pdf=PdfReader(io.BytesIO(raw),strict=False)
    if len(pdf.pages)!=21:
        raise ValueError("original supplemental page count drift")
    page=pdf.pages[TABLE_S1_PAGE_1BASED-1]
    layout=page.extract_text(extraction_mode="layout") or ""
    report=inspect_s1_layout(layout)
    print(json.dumps({
        "status":"PENSTEMON_2026_ORIGINAL_TABLE_S1_HEADER_ONLY",
        "pdf_sha256":ORIGINAL_PDF_SHA,
        "original_page_1based":TABLE_S1_PAGE_1BASED,
        "page_layout_metadata":report,
        "no_SAMN_to_author_voucher_join":True,
    },indent=2,sort_keys=True,ensure_ascii=False))


if __name__=="__main__":
    main()
