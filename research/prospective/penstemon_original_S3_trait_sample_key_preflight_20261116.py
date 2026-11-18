"""Penstemon original Supplement Table S3 sample-ID schema check, no trait values.

The 2026 original PDF SHA256 and 21-page provenance are pinned.
Read only headers and original specimen-key strings from Table S3 page 13;
do not decode or report phenotype, reflectance, nectar or genotype values.
"""
from __future__ import annotations
import hashlib
import io
import json
import re
import sys
from pathlib import Path
import fitz
from pypdf import PdfReader

HERE=Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0,str(HERE))
from penstemon_original_europepmc_table_s1_index_20261115 import source_get, extract_original_pdf_member
from penstemon_original_full52_table_s1_phenotyped_gate_20261115 import PDF_SHA256

ID=re.compile(r"(?<!\w)([A-Z]{2,4})\s*-\s*(\d{2,6})(?:\*)?\b",re.I)
HEAD=re.compile(r"Table\s+S3",re.I)
END=re.compile(r"Table\s+S4",re.I)

def page_metadata(text:str, kind:str)->dict:
    m=HEAD.search(text)
    if not m:raise ValueError(f"Original {kind} page13 Table S3 heading missing")
    sub=text[m.end():]
    end=END.search(sub)
    if end:sub=sub[:end.start()]
    lines=[l.strip() for l in sub.splitlines() if l.strip()]
    named=sorted({f"{x.group(1).upper()}-{x.group(2)}" for x in ID.finditer(sub)})
    # Only short metadata header fragments, with all numerals replaced, are
    # emitted. No floating-point or other quantitative trait fields printed.
    preview=[]
    for l in lines[:14]:
        if not ID.search(l):
            words=re.sub(r"(?<![A-Za-z])\d+(?:\.\d+)?","[number]",l)
            preview.append(words[:130])
    return {
      "extraction_engine":kind,
      "nonblank_text_line_count":len(lines),
      "distinct_voucher_tokens":len(named),
      "voucher_tokens":named,
      "leading_text_sanitized_metadata":preview,
      "no_numeric_trait_values_exported":True,
    }

def main():
    raw=extract_original_pdf_member(source_get())
    if hashlib.sha256(raw).hexdigest()!=PDF_SHA256:
        raise ValueError("Original authoritative Penstemon supplement PDF SHA mismatch")
    a=PdfReader(io.BytesIO(raw),strict=False)
    b=fitz.open(stream=raw,filetype="pdf")
    if len(a.pages)!=21 or len(b)!=21:
        raise ValueError("Original paper supplemental PDF page count drift")
    result={
      "pypdf":page_metadata(a.pages[12].extract_text() or "","pypdf"),
      "pymupdf_sorted":page_metadata(b[12].get_text("text",sort=True),"pymupdf_sorted"),
    }
    print(json.dumps({
      "status":"ORIGINAL_PENSTEMON_TABLE_S3_SOURCE_KEY_SCHEMA_ONLY",
      "original_pdf_sha256":PDF_SHA256,
      "original_page_1based":13,
      "reported_phenotyped_accessions":34,
      "source_engines":result,
      "source_TableS1_52_34_discrepancy_remains_open":True,
      "no_flower_trait_values_or_fitness_genotypes_opened":True,
    },indent=2,sort_keys=True,ensure_ascii=False))

if __name__=="__main__":
    main()
