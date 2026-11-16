"""Render only authentic 2026 Penstemon Table S1 supplemental page for visual QC.

Source PDF sha256 locked and fetched solely from original EuropePMC archive.
Output is temporary run artifact for manual examination, NOT a reconstructed
52-voucher table and NOT a prediction result. No floral numeric data parsed.
"""
from __future__ import annotations
import hashlib
import io
import json
from pathlib import Path
import sys

import fitz

HERE=Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0,str(HERE))

from penstemon_original_europepmc_table_s1_index_20261115 import source_get,extract_original_pdf_member

SHA="3682720dad4f0621bebef3be7014aa667067c5f0dce9160ca8f4c4a758e7f943"
OUT=Path("/tmp/penstemon_2026_S1_source_visual")
def main():
    original=extract_original_pdf_member(source_get())
    if hashlib.sha256(original).hexdigest()!=SHA:
        raise ValueError("Original EuropePMC supplemental PDF checksum mismatch")
    doc=fitz.open(stream=original,filetype="pdf")
    if len(doc)!=21:raise ValueError("Original supplemental source 21 pages expected")
    OUT.mkdir(parents=True,exist_ok=True)
    renders=[]
    # Authoritative Table S1 page and very top of Table S2 page; the latter
    # diagnoses whether a continued S1 row precedes the S2 heading.
    areas=[
      (11,fitz.Rect(doc[10].rect),"Table_S1_page11"),
      (12,fitz.Rect(0,0,doc[11].rect.width,doc[11].rect.height*0.26),"Page12_pre_S2_top"),
    ]
    for page,rect,label in areas:
        pix=doc[page-1].get_pixmap(matrix=fitz.Matrix(2,2),clip=rect,colorspace=fitz.csGRAY,alpha=False)
        target=OUT/(label+".png")
        pix.save(target)
        renders.append({"pdf_page_1based":page,"label":label,
            "render_pixel_width":pix.width,"render_pixel_height":pix.height,
            "png_sha256":hashlib.sha256(target.read_bytes()).hexdigest(),
            "original_pdf_sha256":SHA})
    print(json.dumps({"status":"VERIFIED_ORIGINAL_S1_TABLE_VISUAL_RENDERS_ONLY",
                      "images":renders,
                      "original_pdf_responses_not_fitted":True,
                      "do_not_rescue_52_without_explicit_source_proof":True},
                     indent=2,sort_keys=True))

if __name__=="__main__":main()
