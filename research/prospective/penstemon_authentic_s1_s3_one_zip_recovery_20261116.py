"""Recover exact original Penstemon 2026 Table S1/S3 source layout once.

The official Europe PMC source sometimes stalls while waiting for the final
HTTP chunk. Read bounded ZIP chunks and stop as soon as a complete authentic
ZIP yields the previously certified original PDF SHA256. Never alter source
bytes, silently impute the missing 52nd Table S1 record, or open trait values.
Renders original 11th page and top 12th page for manual source inspection.
"""
from __future__ import annotations
import collections
import hashlib
import io
import json
import re
import sys
import urllib.request
import urllib.error
import zipfile
from pathlib import Path

import fitz
from pypdf import PdfReader

HERE=Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0,str(HERE))
from penstemon_original_europepmc_table_s1_index_20261115 import (
    URL,SOURCE_ZIP_MAX,extract_original_pdf_member
)
from penstemon_original_full52_table_s1_phenotyped_gate_20261115 import (
    PDF_SHA256,ROW,NCBI_source_taxa,identify_original_source_species,origin_phenotyped_geometry
)
from penstemon_original_S3_trait_sample_key_preflight_20261116 import page_metadata

OUT=Path("/tmp/penstemon_authentic_source_recovery_20261116")
CHUNK_SIZE=65_536
MAX_ATTEMPTS=2

def verified_pdf_from_stream():
    attempts=[]
    for attempt in range(1,MAX_ATTEMPTS+1):
        buf=bytearray()
        req=urllib.request.Request(URL,headers={
            "User-Agent":"TTF-Penstemon-original-PDF-recovery/1.0",
            "Accept":"application/zip",
        })
        try:
            with urllib.request.urlopen(req,timeout=85) as response:
                while len(buf)<=SOURCE_ZIP_MAX:
                    part=response.read(min(CHUNK_SIZE,SOURCE_ZIP_MAX+1-len(buf)))
                    if not part:break
                    buf.extend(part)
                    # A valid file must pass original ZIP members, original PDF
                    # file-size AND stable original PDF SHA256. Ignore apparent
                    # ZIP signatures in compressed payloads until those checks.
                    if b"PK\x05\x06" in buf[-(CHUNK_SIZE+128):]:
                        try:pdf=extract_original_pdf_member(bytes(buf))
                        except (zipfile.BadZipFile, ValueError, KeyError, EOFError, OSError):
                            # Source may still be streaming; never accept
                            # partially written ZIP or an unauthenticated PDF.
                            continue
                        if hashlib.sha256(pdf).hexdigest()!=PDF_SHA256:
                            raise ValueError("source original published PDF SHA drift")
                        return pdf,attempts+[{"attempt":attempt,"zip_received":len(buf),"stop":"original-PDF-SHA-success"}]
                # Source may close before an EOCD detection; validate once.
                pdf=extract_original_pdf_member(bytes(buf))
                if hashlib.sha256(pdf).hexdigest()!=PDF_SHA256:
                    raise ValueError("original PDF hash mismatch")
                return pdf,attempts+[{"attempt":attempt,"zip_received":len(buf),"stop":"EOF-and-sha-success"}]
        except (TimeoutError,urllib.error.URLError) as e:
            attempts.append({"attempt":attempt,"status":"SOURCE_STREAM_FAILED","error_type":type(e).__name__,
                             "zip_received_before_stop":len(buf)})
        except ValueError:
            raise
    return None,attempts

def original_s1_partial_source(text:str,ncbi:list[str]):
    rows=[]
    for line in text.splitlines():
        m=ROW.match(line.strip())
        if not m:continue
        hits=identify_original_source_species(m.group("rest"),ncbi)
        rows.append({
          "original_voucher_population_id":m.group("key").upper(),
          "source_species_taxon_match":hits[0] if len(hits)==1 else None,
          "phenotyped":m.group("phenotyped").lower()=="yes",
        })
    ids=[r["original_voucher_population_id"] for r in rows]
    return {
      "status":"INCOMPLETE_SOURCE_DO_NOT_REDEFINE_52_34_GATE",
      "complete_rows_from_pdf":len(rows),
      "unique_source_row_ids":len(set(ids)),
      "explicit_phenotyped_yes":sum(r["phenotyped"] for r in rows),
      "source_taxon_unmatched_ids":[r["original_voucher_population_id"] for r in rows if not r["source_species_taxon_match"]],
      "four_origin_partial_lower_bound_coverage":origin_phenotyped_geometry(rows),
    }

def main():
    OUT.mkdir(parents=True,exist_ok=True)
    pdf,attempts=verified_pdf_from_stream()
    if pdf is None:
        status={
          "status":"ORIGINAL_SUPPLEMENT_TEMPORARILY_UNAVAILABLE",
          "source_url":URL,"failed_attempts":attempts,
          "authentic_pdf_required_sha256":PDF_SHA256,
          "trait_outcomes_not_read":True,
          "not_a_source_missing_voucher_conclusion":True,
        }
        (OUT/"source_recovery_stopped.json").write_text(json.dumps(status,indent=2,sort_keys=True)+"\n")
        print(json.dumps(status,indent=2,sort_keys=True))
        return
    a=PdfReader(io.BytesIO(pdf),strict=False)
    b=fitz.open(stream=pdf,filetype="pdf")
    if len(a.pages)!=21 or len(b)!=21:raise ValueError("unexpected PDF page count")
    page11=a.pages[10].extract_text() or ""
    page12=a.pages[11].extract_text() or ""
    cut=re.search(r"\bTable\s*S2\b",page12,re.I)
    if cut is None:raise ValueError("source page12 Table S2 heading missing")
    original_s1=original_s1_partial_source(page11+"\n"+page12[:cut.start()],sorted(NCBI_source_taxa()))
    s3=page_metadata(a.pages[12].extract_text() or "","pypdf")
    renders=[]
    for pageno,clip,name in [
        (11,b[10].rect,"original_table_S1_page11"),
        (12,fitz.Rect(0,0,b[11].rect.width,b[11].rect.height*.3),"original_page12_top_before_S2"),
    ]:
        pix=b[pageno-1].get_pixmap(matrix=fitz.Matrix(1.9,1.9),
            colorspace=fitz.csGRAY,alpha=False,clip=clip)
        fn=OUT/(name+".png")
        pix.save(fn)
        renders.append({"source_pdf_page_1based":pageno,"image_file_name":fn.name,
             "image_sha256":hashlib.sha256(fn.read_bytes()).hexdigest(),
             "image_pixels":[pix.width,pix.height]})
    result={
      "status":"AUTHENTIC_PENSTEMON_S1_S3_SOURCE_VISUAL_AND_MISSINGNESS_DIAGNOSIS",
      "original_pdf_sha256":PDF_SHA256,
      "network_source_attempts":attempts,
      "S1":original_s1,
      "S3_identity_schema":s3,
      "original_source_renders":renders,
      "not_52_genome_34_flower_link_complete":True,
      "no_quantitative_floral_or_genomic_effects_opened":True,
    }
    receipt=OUT/"original_source_audit.json"
    receipt.write_text(json.dumps(result,indent=2,sort_keys=True,ensure_ascii=False)+"\n")
    # Print aggregate only, not original pages or trait values.
    print(json.dumps({
      "status":result["status"],
      "S1_complete_rows":original_s1["complete_rows_from_pdf"],
      "S1_source_yes_flags":original_s1["explicit_phenotyped_yes"],
      "S1_taxon_unmatched_ids":original_s1["source_taxon_unmatched_ids"],
      "S1_four_origin_partial_coverage":original_s1["four_origin_partial_lower_bound_coverage"],
      "S3_source_voucher_id_tokens":s3["distinct_voucher_tokens"],
      "S3_leading_sanitized_metadata":s3["leading_text_sanitized_metadata"][:8],
      "rendered_images":renders,"no_phenotype_effect_values_read":True,
    },indent=2,sort_keys=True,ensure_ascii=False))

if __name__=="__main__":main()
