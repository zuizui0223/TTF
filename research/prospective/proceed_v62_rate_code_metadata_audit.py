"""Inspect exact published PROCEED v6.2 BASIC SOURCE CODE only, never run it.

Fetches public R/text source from Borealis Dataverse and verifies official
file ID and MD5. No phenotype, fitness, genotype or rate-data record is read.
Selected short source-code lines are diagnostic provenance, not a modeled
outcome. No reuse or reopening of frozen TTF response families.
"""
from __future__ import annotations
import argparse
import hashlib
import json
import re
import urllib.request
from pathlib import Path

from proceed_v62_borealis_header_only import (
    validated_get_catalog, verify_frozen_catalog
)
HERE=Path(__file__).resolve().parent
CONTRACT=HERE/"proceed_v62_rate_code_metadata_contract_v01.json"
SOURCE_MANIFEST=HERE/"proceed_v62_source_header_contract_v01.json"

def code_signals(body,terms,max_lines=35):
    # Only read public code comments/expressions, not source data tables.
    txt=body.decode("utf-8-sig","replace")
    lines=txt.splitlines()
    matched=[]
    scores={}
    for i,line in enumerate(lines,1):
        stripped=line.strip()
        if not stripped or stripped.startswith("#") and len(stripped)<5:
            continue
        labels=[term for term in terms if term.lower() in stripped.lower()]
        if not labels:
            continue
        score=sum(4 if term in ("haldane","darwin","log10","gtime","generations") else 1
                  for term in labels)
        score+=int(bool(re.search(r"=|<-|function|/",stripped)))*2
        scores[i]=score
        matched.append({"line_number":i,"source_text_excerpt":stripped[:170],
                        "controlled_terms":labels,"score":score})
    # Key lines in source order selected by predeclared static term search;
    # no source response rates or performance selected these.
    preferred=sorted(matched,key=lambda x:(-x["score"],x["line_number"]))[:max_lines]
    preferred.sort(key=lambda x:x["line_number"])
    # Explicit codebook arithmetic snippets only, not a data table. Limited
    # to a bounded original BasicCode section describing Darwin/Haldane
    # admissibility, numerator/denominator and published transforms.
    targeted=[]
    for i,line in enumerate(lines,1):
        if not (20<=i<=132):continue
        p=line.strip()
        if (p.startswith("#") and not any(k in p.lower() for k in
              ("darwin","haldane","generation","transf_data","pooled"))) or not p:
            continue
        if (any(token in p for token in ("datos$D[","datos$H[","datos$numD[",
            "datos$numH[","datos$a_D","datos$a_H","datos$pSD[",
            "datos$La_D","datos$La_H","datos$okD[","datos$okH[",
            "okD_ln","okH","D <-","H <-")) or
            (p.startswith("#") and any(k in p.lower() for k in
             ("darwin","haldane","generation","transf_data","pooled")))):
            targeted.append({"line":i,"code_excerpt":p[:180]})
    targeted=targeted[:65]
    return {
        "line_count":len(lines),
        "bounded_rate_unit_arithmetic_code_excerpts":targeted,
        "number_of_source_code_lines_matching_declared_terms":len(matched),
        "first_controlled_term_matches":preferred,
        "whole_source_code_exported":False,
        "code_executed":False,
        "genetic_trait_rate_rows_read":False,
        "bio_response_fit_or_studied":False
    }


def run(args):
    contract=json.loads(CONTRACT.read_text())
    if contract.get("status")!="PRE_RESPONSE_PUBLIC_ANALYSIS_CODE_ONLY_AUDIT_FROZEN":
        raise ValueError("Public code-only probe not frozen")
    source=contract["original_code_file"]
    frozen=json.loads(SOURCE_MANIFEST.read_text())["current_source"]
    out={"schema":"ttf_proceed_v62_exact_source_code_unit_definition_audit_v01",
         "source_id":source["id"],"source_name":source["name"],
         "genetic_trait_rate_rows_read":False,
         "bio_response_fit_or_studied":False,
         "code_executed":False}
    try:
        catalog=validated_get_catalog(
            "https://borealisdata.ca/api/datasets/:persistentId/?persistentId=doi%3A10.5683%2FSP3%2FNXSL3Q"
        )
        verify_frozen_catalog(catalog,frozen)
        entries=catalog["data"]["latestVersion"]["files"]
        meta=[z for z in entries if z.get("dataFile",{}).get("id")==source["id"]]
        if len(meta)!=1:
            raise ValueError("Original published BasicCode ID absent/duplicated")
        df=meta[0]["dataFile"]
        if (df.get("filename")!=source["name"] or
            df.get("md5")!=source["md5"] or
            df.get("filesize")!=source["catalog_size_bytes"] or
            meta[0].get("restricted")):
            raise ValueError("Original BasicCode source metadata MD5/size/name changed")
        url=source["source_endpoint"]
        if url!="https://borealisdata.ca/api/access/datafile/1044644?format=original":
            raise ValueError("Public code source endpoint drift")
        req=urllib.request.Request(url,headers={
            "User-Agent":"TTF-PROCEED-source-code-read-only-v01"
        })
        with urllib.request.urlopen(req,timeout=60) as response:
            data=response.read(100000)
        md5=hashlib.md5(data).hexdigest()
        sha=hashlib.sha256(data).hexdigest()
        if md5!=source["md5"]:
            raise ValueError("Original BasicCode native MD5 not reproducible; STOP")
        report=code_signals(data,contract["controlled_terms"])
        out.update(report)
        out["original_file_md5_verified"]=md5
        out["original_file_sha256"]=sha
        out["original_file_bytes"]=len(data)
        out["catalog_ingested_size_bytes"]=source["catalog_size_bytes"]
        out["status"]="ORIGINAL_CODE_IDENTITY_VERIFIED_SHORT_FORMULA_EXCERPTS_ONLY"
    except Exception as exc:
        out["status"]="SOURCE_CODE_UNAVAILABLE_NOT_BIOLOGICAL_RESULT"
        out["error_type"]=type(exc).__name__
        out["error_message"]=str(exc)[:250]
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(out,sort_keys=True,indent=2)+"\n")
    print(json.dumps(out,sort_keys=True,ensure_ascii=False))


if __name__=="__main__":
    ap=argparse.ArgumentParser()
    ap.add_argument("--output",required=True,type=Path)
    run(ap.parse_args())
