"""Diagnose missing 52nd Penstemon original Table S1 voucher visually.

Original authenticated Europe PMC 2026 supplement PDF sha256 fixed.
Render ONLY PDF page 11 (Table S1) to PNG for manual original-source
verification, plus a bounded identity-only text-line diagnostic.
No Table S2/S3 values, no floral quantitative outcomes or genotypes.
"""
from __future__ import annotations
import hashlib
import io
import json
import os
import re
from pathlib import Path

import fitz
from pypdf import PdfReader
import sys

HERE=Path(__file__).resolve().parent
if str(HERE) not in sys.path:sys.path.insert(0,str(HERE))
from penstemon_original_europepmc_table_s1_index_20261115 import (
    source_get,extract_original_pdf_member
)
from penstemon_original_full52_table_s1_phenotyped_gate_20261115 import (
    PDF_SHA256,ROW,MAYBE_ID
)
PAGE=11
OUTPUT=Path("/tmp/penstemon_2026_original_table_s1_page11.png")


def source_line_audit(text:str)->dict:
    if "Table S1" not in text:
        raise ValueError("original source title absent")
    lines=[x.strip() for x in text.splitlines() if x.strip()]
    audited=[]
    for i,raw in enumerate(lines,1):
        is_full=bool(ROW.match(raw))
        looks=bool(MAYBE_ID.match(raw))
        ends=bool(re.search(r"\b(?:yes|no)\s*$",raw,re.I))
        audited.append({
           "line_number":i,
           "classified_as_complete_voucher":is_full,
           "looks_like_original_voucher_id":looks,
           "ends_in_phenotyped_yes_no":ends,
           "original_metadata_line_preview":raw[:165],
        })
    return {
       "nonblank_source_page11_line_count":len(lines),
       "complete_voucher_lines":sum(x["classified_as_complete_voucher"] for x in audited),
       "unmatched_lines_ending_yes_no":[x for x in audited if x["ends_in_phenotyped_yes_no"] and not x["classified_as_complete_voucher"]],
       "lines_look_like_voucher_but_not_complete":[x for x in audited if x["looks_like_original_voucher_id"] and not x["classified_as_complete_voucher"]],
       "source_lines_preview_for_visual_diagnosis":audited,
       "not_a_phenotype_response_analysis":True,
    }


def main():
    raw=extract_original_pdf_member(source_get())
    if hashlib.sha256(raw).hexdigest()!=PDF_SHA256:
        raise ValueError("source original PDF SHA256 drift")
    doc=PdfReader(io.BytesIO(raw),strict=False)
    if len(doc.pages)!=21:raise ValueError("source original PDF page count drift")
    text=doc.pages[PAGE-1].extract_text() or ""
    audit=source_line_audit(text)
    image=fitz.open(stream=raw,filetype="pdf")
    if len(image)!=21:raise ValueError("PyMuPDF source page count inconsistent")
    pix=image[PAGE-1].get_pixmap(matrix=fitz.Matrix(1.5,1.5),alpha=False)
    pix.save(str(OUTPUT))
    pdf_digest=hashlib.sha256(raw).hexdigest()
    print(json.dumps({
       "status":"ORIGINAL_TABLE_S1_VISUAL_DIAGNOSIS_SOURCE_ONLY",
       "pdf_sha256":pdf_digest,
       "page_1based":PAGE,
       "render_png_path":str(OUTPUT),
       "render_png_sha256":hashlib.sha256(OUTPUT.read_bytes()).hexdigest(),
       "render_pixels":[pix.width,pix.height],
       "line_audit":audit,
       "no_phenotypic_magnitude_or_genotypes_read":True,
       "requires_original_image_review_before_rescuing_52nd_voucher":True,
    },indent=2,sort_keys=True,ensure_ascii=False))


if __name__=="__main__":main()
