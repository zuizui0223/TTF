"""Penstemon 2026 ORIGINAL Table S1 52 voucher IDs × phenotyped flags only.

The source-protocol was committed before this original all-row presence audit:
PENSTEMON_FULL_S1_PHENOTYPED_VOUCHER_SOURCE_RULE_20261115.md

PDF: Europe PMC PMC13419788, 21 pages, page 11 Table S1,
SHA256 3682720dad4f0621bebef3be7014aa667067c5f0dce9160ca8f4c4a758e7f943.
No floral quantitative effects, spectra, genotypes or evolutionary scores read.
"""
from __future__ import annotations
import collections
import hashlib
import io
import json
import re
import sys
from pathlib import Path
from pypdf import PdfReader

HERE=Path(__file__).resolve().parent
if str(HERE) not in sys.path:sys.path.insert(0,str(HERE))
from penstemon_original_europepmc_table_s1_index_20261115 import (
    source_get,extract_original_pdf_member
)
from penstemon_four_origin_source_genome_blocks_20261115 import BLOCKS

PDF_SHA256="3682720dad4f0621bebef3be7014aa667067c5f0dce9160ca8f4c4a758e7f943"
TABLE_PAGE_1BASED=11
GENOME_SOURCE=HERE/"penstemon_habroanthus_52_ncbi_biosample_taxon_receipt_20261115.json"
GENOME_SOURCE_BLOB="acdf2997e4d338fcb879d0323dc4c6a54f4fe413"
# This is an original source *voucher ID* prefix. A SAMN is NOT this ID.
ROW=re.compile(r"^(?P<key>[A-Z]{2,4}-\d{2,6})(?P<star>\*)?\s+(?P<rest>.*\S)\s+(?P<phenotyped>yes|no)\s*$",re.I)
MAYBE_ID=re.compile(r"^[A-Z]{2,4}-\d{2,6}\*?\b",re.I)


def git_blob_sha(raw:bytes)->str:
    return hashlib.sha1(f"blob {len(raw)}\0".encode()+raw).hexdigest()


def NCBI_source_taxa() -> dict[str,int]:
    raw=GENOME_SOURCE.read_bytes()
    if git_blob_sha(raw)!=GENOME_SOURCE_BLOB:
        raise ValueError("original 52-BioSample source receipt Git blob drift")
    d=json.loads(raw)
    s=d.get("source_matched_taxonomy") or {}
    if s.get("status")!="ALL_52_ORIGINAL_BIOSAMPLES_TAXON_SOURCE_MATCHED":
        raise ValueError("52 original NCBI samples not authenticated")
    counts=s.get("source_taxon_sample_counts") or {}
    if sum(counts.values())!=52 or len(counts)!=21:
        raise ValueError("original NCBI source taxon coverage changed")
    return counts


def identify_original_source_species(line:str,ncbi_species:list[str]):
    # The PDF extracts intraword whitespace in scientific names, e.g.,
    # P. speci osus. Only whitespace, not synonymy or spelling, is collapsed.
    compact=re.sub(r"\s+","",line.lower())
    match=[]
    for full in ncbi_species:
        if not full.startswith("Penstemon "):continue
        needle="p."+full[len("Penstemon "):].lower().replace(" ","")
        if needle in compact:
            match.append(full)
    return match


def parse_s1_page(text:str,ncbi_species:list[str],expected_n:int=52,expected_yes:int=34)->dict:
    if not re.search(r"\bTable\s*S1\b",text,re.I):
        raise ValueError("Table S1 heading not in original source page")
    if not re.search(r"\bPhenotyped\b",text,re.I):
        raise ValueError("Original table phenotyped column not present")
    rows=[]
    suspicious=[]
    for line in text.splitlines():
        cell=line.strip()
        if not cell:continue
        if not MAYBE_ID.match(cell):
            continue
        m=ROW.match(cell)
        if not m:
            suspicious.append(cell[:180])
            continue
        possible=identify_original_source_species(m.group("rest"),ncbi_species)
        rows.append({
            "original_voucher_population_id":m.group("key").upper(),
            "original_5taxon_introgression_marker_asterisk":bool(m.group("star")),
            "source_species_taxon_match":possible[0] if len(possible)==1 else None,
            "species_match_candidate_count":len(possible),
            "phenotyped":m.group("phenotyped").lower()=="yes",
            # No original collection location or numerical outcome exported.
        })
    ids=[x["original_voucher_population_id"] for x in rows]
    n_yes=sum(x["phenotyped"] for x in rows)
    if len(rows)!=expected_n or len(set(ids))!=expected_n or suspicious:
        raise ValueError(
          f"Original Table S1 row structure not fully resolved: accepted={len(rows)} "
          f"unique={len(set(ids))} suspicious={len(suspicious)} "
          f"first_unparsed={suspicious[:3]}"
        )
    if n_yes!=expected_yes:
        raise ValueError(f"Original source phenotyped count {n_yes} differs from paper {expected_yes}")
    species_counts=dict(sorted(collections.Counter(
        x["source_species_taxon_match"] for x in rows if x["source_species_taxon_match"]
    ).items()))
    sample_status=collections.Counter(
        "phenotyped" if x["phenotyped"] else "not_phenotyped" for x in rows
    )
    unmatched=[x["original_voucher_population_id"] for x in rows if not x["source_species_taxon_match"]]
    return {
        "status":"ORIGINAL_PDF_TABLE_S1_52_AUTHOR_VOUCHERS_34_PHENOTYPED",
        "source_sample_count":len(rows),
        "source_phenotyped_count":n_yes,
        "original_voucher_rows":rows,
        "source_taxon_matched_count":len(rows)-len(unmatched),
        "unmatched_original_voucher_ids_no_taxon_synonym_guess":unmatched,
        "source_species_match_counts":species_counts,
        "no_quantitative_floral_or_genetic_values_read":True,
        "sample_to_SAMN_individual_ID_join":"NOT_VERIFIED",
        "source_phenotyped_flag_is_not_a_trait_value":True,
    }


def origin_phenotyped_geometry(rows:list[dict])->dict:
    results={}
    for origin,group in BLOCKS.items():
        species=[group["derived"]]+list(group["sisters"])
        srows=[r for r in rows if r["source_species_taxon_match"] in species]
        counts={t:{
           "original_source_voucher_count":sum(r["source_species_taxon_match"]==t for r in srows),
           "original_source_phenotyped_count":sum(
               r["source_species_taxon_match"]==t and r["phenotyped"] for r in srows
           )
        } for t in species}
        results[origin]={
           "source_species_phenotyped_by_original_taxon":counts,
           "source_block_voucher_count":len(srows),
           "source_block_phenotyped_count":sum(r["phenotyped"] for r in srows),
           "origin_both_derived_and_any_bee_phenotyped":(
               counts[group["derived"]]["original_source_phenotyped_count"]>0
               and any(counts[s]["original_source_phenotyped_count"]>0 for s in group["sisters"])
           ),
           "same_NCBI_SAMN_individuals_verified":False,
        }
    return results


def main():
    data=extract_original_pdf_member(source_get())
    if hashlib.sha256(data).hexdigest()!=PDF_SHA256:
        raise ValueError("original PDF page source SHA256 drift")
    pdf=PdfReader(io.BytesIO(data),strict=False)
    if len(pdf.pages)!=21:raise ValueError("original supplemental PDF page count drift")
    text=pdf.pages[TABLE_PAGE_1BASED-1].extract_text() or ""
    ncbi=NCBI_source_taxa()
    page=parse_s1_page(text,sorted(ncbi))
    counts=page["source_species_match_counts"]
    mismatch={
       t:{"n_original_2026_table_s1_vouchers":counts.get(t,0),"n_original_2026_ncbi_SAMN":ncbi[t]}
       for t in sorted(ncbi) if counts.get(t,0)!=ncbi[t]
    }
    blocks=origin_phenotyped_geometry(page["original_voucher_rows"])
    print(json.dumps({
        "status":"PENSTEMON_ORIGINAL_2026_S1_PUBLISHED_PHENOTYPED_KEYS_ONLY",
        "original_paper_doi":"10.64898/2026.07.21.739878",
        "original_pdf_sha256":PDF_SHA256,
        "original_source_table_1based_page":TABLE_PAGE_1BASED,
        "original_voucher_and_phenotyped_source":page,
        "origin_phenotyped_presence":blocks,
        "source_taxon_count_disagreement_NCBI_vs_PDF":mismatch,
        "all_four_origins_both_types_phenotyped_at_paper_species_label_level":
            all(x["origin_both_derived_and_any_bee_phenotyped"] for x in blocks.values()),
        "sample_to_SAMN_individual_join":"NOT_VERIFIED",
        "TTF_P_predictive_effect_computed":False,
        "TTF_G_or_RI_effect_computed":False,
    },indent=2,sort_keys=True,ensure_ascii=False))


if __name__=="__main__":main()
