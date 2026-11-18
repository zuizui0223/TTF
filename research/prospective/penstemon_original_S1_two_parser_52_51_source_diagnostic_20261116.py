"""Independent two-parser identity-only diagnosis of Penstemon original Table S1.

Source image/PDF contents SHA-pinned. Reads only page 11 and page-12 text
before Table S2; never exports floral trait values or genotype outcomes.
This diagnoses whether the published 52 versus pypdf 51 sample discrepancy
is due to PDF text order, line splitting or source table omission.
"""
from __future__ import annotations
import collections
import hashlib
import io
import json
import re
from pathlib import Path
import sys

import fitz
from pypdf import PdfReader

HERE=Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0,str(HERE))

from penstemon_original_europepmc_table_s1_index_20261115 import source_get,extract_original_pdf_member
from penstemon_original_full52_table_s1_phenotyped_gate_20261115 import (
    PDF_SHA256,ROW,MAYBE_ID,identify_original_source_species,NCBI_source_taxa
)

ID_ANYWHERE=re.compile(r"(?<![A-Za-z0-9])(?:(?P<unhyphen>TTC026503)|(?P<prefix>[A-Z]{2,4})\s*-\s*(?P<digits>\d{2,6})(?P<star>\*)?)",re.I)
TABLE_S2=re.compile(r"\bTable\s*S2\b",re.I)


def clean_id(match):
    if match.group("unhyphen"):
        return "TTC026503"  # exact visually verified source voucher, not a general ID grammar
    return match.group("prefix").upper()+"-"+match.group("digits")


def text_profile(text:str,source:str,ncbi_species:list[str])->dict:
    lines=[line.strip() for line in text.splitlines() if line.strip()]
    exact={}
    id_mentions=collections.Counter()
    missing_candidate=[]
    phenotyped_yes=0
    phenotyped_no=0
    species_count=collections.Counter()
    for line in lines:
        for m in ID_ANYWHERE.finditer(line):
            id_mentions[clean_id(m)]+=1
        source_match=ROW.match(line)
        if source_match:
            key=source_match.group("key").upper()
            if key in exact:
                raise ValueError(f"duplicate complete published voucher row: {key}")
            available=source_match.group("phenotyped").lower()=="yes"
            phenotyped_yes+=int(available)
            phenotyped_no+=int(not available)
            species=identify_original_source_species(source_match.group("rest"),ncbi_species)
            exact[key]={"phenotyped":available,"source_taxon":species[0] if len(species)==1 else None,
                        "n_source_taxon_candidates":len(species)}
            if len(species)==1:species_count[species[0]]+=1
        elif MAYBE_ID.match(line):
            v=ID_ANYWHERE.search(line)
            missing_candidate.append(clean_id(v) if v else "unparseable")
    return {
      "source":source,
      "nonblank_lines":len(lines),
      "exact_rows":len(exact),
      "exact_phenotyped_yes":phenotyped_yes,
      "exact_phenotyped_no":phenotyped_no,
      "unique_identifiers_any_position":len(id_mentions),
      "id_mentions_not_complete_rows":sorted(set(id_mentions)-set(exact)),
      "complete_rows_not_found_by_anywhere_regex":sorted(set(exact)-set(id_mentions)),
      "incomplete_original_id_prefixed_lines":missing_candidate,
      "exact_row_ids":sorted(exact),
      "source_taxa_from_complete_rows":dict(sorted(species_count.items())),
      "unmatched_source_taxon_ids":sorted(k for k,v in exact.items() if not v["source_taxon"]),
      "no_trait_magnitudes_read_or_exported":True,
    }


def select_s1_source_text(page11:str,page12:str)->str:
    if not re.search(r"\bTable\s*S1\b",page11,re.I):
        raise ValueError("original supplemental page 11 Table S1 absent")
    split=TABLE_S2.search(page12)
    if not split:raise ValueError("original supplemental page 12 Table S2 absent")
    return page11+"\n"+page12[:split.start()]


def main():
    raw=extract_original_pdf_member(source_get())
    if hashlib.sha256(raw).hexdigest()!=PDF_SHA256:
        raise ValueError("original 2026 Penstemon supplement PDF SHA256 drift")
    a=PdfReader(io.BytesIO(raw),strict=False)
    b=fitz.open(stream=raw,filetype="pdf")
    if len(a.pages)!=21 or len(b)!=21:raise ValueError("original page count mismatch")
    names=sorted(NCBI_source_taxa())
    versions={
      "pypdf_plain":select_s1_source_text(a.pages[10].extract_text() or "",a.pages[11].extract_text() or ""),
      "pypdf_layout":select_s1_source_text(
         (a.pages[10].extract_text(extraction_mode="layout") or "") if "Table S1" in (a.pages[10].extract_text(extraction_mode="layout") or "") else (a.pages[10].extract_text() or ""),
         a.pages[11].extract_text() or ""),
      "pymupdf_plain":select_s1_source_text(b[10].get_text("text",sort=False),b[11].get_text("text",sort=False)),
      "pymupdf_sorted":select_s1_source_text(b[10].get_text("text",sort=True),b[11].get_text("text",sort=True)),
    }
    results={key:text_profile(text,key,names) for key,text in versions.items()}
    p=results["pypdf_plain"]
    compare={
      k:{
        "identifiers_absent_from_pypdf_plain":sorted(set(v["exact_row_ids"])-set(p["exact_row_ids"])),
        "pypdf_plain_identifiers_not_in_this_parser":sorted(set(p["exact_row_ids"])-set(v["exact_row_ids"])),
        "all_position_identifiers_absent_from_pypdf_plain_full_rows":sorted(set(v["exact_row_ids"])-set(p["id_mentions_not_complete_rows"])-set(p["exact_row_ids"])),
      } for k,v in results.items() if k!="pypdf_plain"
    }
    ncbi=NCBI_source_taxa()
    difference={t:{"n_biosample":ncbi[t],"n_pdf_pypdf_complete":p["source_taxa_from_complete_rows"].get(t,0)}
       for t in names if ncbi[t]!=p["source_taxa_from_complete_rows"].get(t,0)}
    print(json.dumps({
      "status":"PENSTEMON_SOURCE_ONLY_52_VS_51_EXTRACTOR_DIAGNOSIS",
      "original_pdf_sha256":PDF_SHA256,
      "published_52_sequenced":52,
      "published_34_phenotyped":34,
      "source_pages_1based":[11,12],
      "per_engine":results,
      "differences_from_pypdf_plain":compare,
      "NCBI_taxa_count_mismatch_with_51_rows":difference,
      "no_experimental_or_genetic_numeric_data_read":True,
      "no_P_TTF_effect_opened":True,
      "do_not_admit_phenotype_sample_join_until_52_34_confirmed":True
    },indent=2,sort_keys=True,ensure_ascii=False))

if __name__=="__main__":
    main()
