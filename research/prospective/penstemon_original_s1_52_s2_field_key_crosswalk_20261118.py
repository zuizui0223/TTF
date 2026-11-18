"""Penstemon 2026 supplement: visually corrected 52 S1 vouchers and original S2 field keys.

STRICT source-only post-visual audit, NOT a publication-faithful 34-sample
measurement join. Original Europe PMC PDF SHA256 is unchanged. Fields
read: original voucher ID, source species string, S1 yes/no, S2 field
voucher ID/species name. No floral trait magnitudes, genotypes or p-values.
"""
from __future__ import annotations

import hashlib
import io
import json
import re
import sys
import time
import urllib.error
from pathlib import Path

from pypdf import PdfReader

HERE=Path(__file__).resolve().parent
if str(HERE) not in sys.path:sys.path.insert(0,str(HERE))

from penstemon_original_europepmc_table_s1_index_20261115 import (
    source_get,extract_original_pdf_member
)
from penstemon_original_full52_table_s1_phenotyped_gate_20261115 import (
    PDF_SHA256, NCBI_source_taxa, parse_s1_page, origin_phenotyped_geometry
)

# S2 field voucher IDs use a different source table from the genomes. Only
# compare normalized exact original alphanumeric collector-voucher token,
# never assume a field specimen is among the 52 sequenced accessions.
S2_ROW=re.compile(r"^\s*(?P<key>[A-Z]{2,4}\s*-?\s*\d{2,6})\s+(?P<species>P\.\s*[A-Za-z]+)\s+(?P<rest>.+)$",re.I)
S2_TITLE=re.compile(r"\bTable\s+S2\b",re.I)
S3_TITLE=re.compile(r"\bTable\s+S3\b",re.I)


def get_original_pdf():
    for attempt in range(3):
        try:
            original=extract_original_pdf_member(source_get())
            if hashlib.sha256(original).hexdigest()!=PDF_SHA256:
                raise ValueError("Original Europe PMC supplemental PDF SHA256 mismatch")
            return original
        except (TimeoutError,urllib.error.URLError):
            if attempt==2:raise
            time.sleep(1+attempt*2)
    raise ValueError("source retrieval failed")


def normalize_voucher(source):
    """Remove separator spaces/hyphens only, never resynthesize numeric IDs."""
    return re.sub(r"[\s-]+","",source).upper()


def parse_s2_field_identity(page_text: str) -> dict:
    match=S2_TITLE.search(page_text)
    if not match:
        raise ValueError("original Table S2 heading absent on page 12")
    trailing=page_text[match.end():]
    end=S3_TITLE.search(trailing)
    if end:trailing=trailing[:end.start()]
    found=[]
    candidate=[]
    for line in trailing.splitlines():
        line=line.strip()
        if not line:continue
        m=S2_ROW.match(line)
        if m:
            key=re.sub(r"\s+","",m.group("key")).upper()
            found.append({
                "original_s2_field_voucher_id":key,
                "species_name_prefix":m.group("species").strip()[:50],
            })
        elif re.match(r"^[A-Z]{2,4}\s*-?\s*\d{2,6}",line,re.I):
            candidate.append(line[:80])
    ids=[r["original_s2_field_voucher_id"] for r in found]
    if len(ids)!=len(set(ids)):
        raise ValueError("duplicate original Table S2 field voucher row")
    if candidate:
        # Fail closed rather than silently assert original S2 is complete.
        return {
           "status":"SOURCE_S2_IDENTITY_PARSER_INCOMPLETE",
           "complete_source_rows":len(found),
           "unparsed_id_prefix_line_count":len(candidate),
           "unparsed_line_first3_redacted":candidate[:3],
           "no_source_numeric_effects_read":True,
        }
    return {
        "status":"SOURCE_S2_FIELD_IDENTITIES_PARSED",
        "source_rows":len(ids),
        "original_field_voucher_ids":ids,
        "original_field_species_name_prefixes":[r["species_name_prefix"] for r in found],
        "no_source_numeric_effects_read":True,
    }


def main():
    pdf=PdfReader(io.BytesIO(get_original_pdf()),strict=False)
    if len(pdf.pages)!=21:raise ValueError("Original PDF 21-page count changed")
    page11=pdf.pages[10].extract_text() or ""
    page12=pdf.pages[11].extract_text() or ""
    s2match=S2_TITLE.search(page12)
    if not s2match:raise ValueError("Table S2 absent at start of original page12")
    # The original Table S1 ends on page11; source no continuation before S2.
    # Recovered visually verified row TTC026503 is exact and Phenotyped=no.
    s1=parse_s1_page(page11,sorted(NCBI_source_taxa()),expected_n=52,expected_yes=30)
    keys={r["original_voucher_population_id"]:r for r in s1["original_voucher_rows"]}
    if "TTC026503" not in keys or keys["TTC026503"]["phenotyped"]:
        raise ValueError("original visually authenticated TTC026503 cardinalis no source row missing")
    if keys["TTC026503"]["source_species_taxon_match"]!="Penstemon cardinalis":
        raise ValueError("original TTC026503 cardinalis name mismatched")
    original_s2=parse_s2_field_identity(page12)
    if original_s2["status"]=="SOURCE_S2_FIELD_IDENTITIES_PARSED":
        normal_s1={normalize_voucher(k):k for k in keys}
        original_s2["s2_field_ids_also_present_as_s1_genome_vouchers"]=[
            {"s2":k,"s1":normal_s1[normalize_voucher(k)]}
            for k in original_s2["original_field_voucher_ids"]
            if normalize_voucher(k) in normal_s1
        ]
        original_s2["s2_field_vouchers_not_in_s1"]=[
            k for k in original_s2["original_field_voucher_ids"]
            if normalize_voucher(k) not in normal_s1
        ]
        original_s2["s2_source_names_are_not_phenotype_values"]=True
    sample_count_s1=s1["source_sample_count"]
    yes_count=s1["source_phenotyped_count"]
    no_count=sample_count_s1-yes_count
    output={
        "status":"S1_VISUAL_RECORD_52_CONFIRMED_S1_YES_30_NOT_34",
        "original_pdf_sha256":PDF_SHA256,
        "visually_authenticated_special_id":"TTC026503",
        "special_id_source_taxon":"Penstemon cardinalis",
        "special_id_phenotyped":False,
        "s1_source_voucher_count":sample_count_s1,
        "s1_source_yes_count":yes_count,
        "s1_source_no_count":no_count,
        "manuscript_reports_flower_phenotyped_accessions":34,
        "original_NCBI_sample_count":sum(NCBI_source_taxa().values()),
        "source_taxon_match_count":s1["source_taxon_matched_count"],
        "unmatched_literal_source_species_keys":s1["unmatched_original_voucher_ids_no_taxon_synonym_guess"],
        "S1_source_yes_by_four_origin":origin_phenotyped_geometry(s1["original_voucher_rows"]),
        "S2_field_source":original_s2,
        "no_S3_trait_values_read":True,
        "no_flower_trait_or_genotype_values_read":True,
        "S1_yes_30_is_not_manuscript_34":True,
        "no_individual_voucher_to_SAMN_mapping_inferred":True,
        "TTF_P_G_RI_effects_not_evaluated":True,
    }
    print(json.dumps(output,indent=2,ensure_ascii=False,sort_keys=True))


if __name__=="__main__":
    main()
