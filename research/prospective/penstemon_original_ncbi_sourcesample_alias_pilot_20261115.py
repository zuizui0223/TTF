"""First-five public Penstemon NCBI source original alias versus 2023 tips.

Fixed metadata-only diagnostic before inspecting 2026 ESummary sourcesample.
No sampledata contents, sequence, genotype, or floral trait values are read.
"""
from __future__ import annotations
import csv
import hashlib
import io
import json
import sys
import urllib.error
import urllib.parse
import urllib.request
from pathlib import Path

HERE=str(Path(__file__).resolve().parent)
if HERE not in sys.path:sys.path.insert(0,HERE)
from penstemon_habroanthus_ncbi_biosample_source_pilot_20261115 import (
    source_samples,numeric_uids_from_esearch,NCBI_SEARCH,NCBI,
    source_metadata_aggregate,MAX_BYTES
)

SOURCE_REPO="lcwheeler/penstemon_novaseq"
SOURCE_COMMIT="fd93699387a8"
SOURCE_FILE="4_analysis/twisst_analysis/target_trios/groups.tsv"
SOURCE_BLOB="360d904346db4ebc17b13d559aed5cd13dde2867"
N=5


def original_2023_mapped_tips()->dict[str,str]:
    url=f"https://raw.githubusercontent.com/{SOURCE_REPO}/{SOURCE_COMMIT}/{SOURCE_FILE}"
    with urllib.request.urlopen(url,timeout=35) as resp:raw=resp.read(15_000)
    sha=hashlib.sha1(f"blob {len(raw)}\0".encode()+raw).hexdigest()
    if sha!=SOURCE_BLOB:
        raise ValueError("Original older public mapping Git blob drift")
    tips={}
    for row in csv.reader(io.StringIO(raw.decode("utf-8-sig")),delimiter="\t"):
        if len(row)!=2 or not all(row) or row[0] in tips:
            raise ValueError("Original 2023 source tip/species map invalid")
        tips[row[0]]=row[1]
    if len(tips)!=18:raise ValueError("2023 original source subset count drift")
    return tips


def inspect_alias_fields(raw:bytes,source_ids:list[str],original_tip_map:dict[str,str])->dict:
    if len(raw)>MAX_BYTES:raise ValueError("NCBI source summary over cap")
    meta=source_metadata_aggregate(raw,source_ids)
    if meta.get("status")!="NCBI_BIOSAMPLE_SOURCE_SUMMARY_ONLY":
        return {"status":"NCBI_ESUMMARY_SOURCE_NOT_EVALUABLE"}
    if len(meta["records"])!=len(source_ids) or set(meta["unmatched_source_accessions"]):
        raise ValueError("Source expected five BioSample IDs not completely matched")
    data=json.loads(raw).get("result") or {}
    rows=[]
    for uid in data.get("uids") or []:
        source=data.get(str(uid)) or {}
        samn=source.get("accession")
        aliases=source.get("sourcesample")
        identifiers=source.get("identifiers")
        alias=aliases.strip() if isinstance(aliases,str) and len(aliases)<257 else None
        matching=alias in original_tip_map if alias else False
        rows.append({
            "biosample_accession":samn,
            "source_submitter_sample_key":alias,
            "original_key_available":bool(alias),
            "original_source_attribute_type":type(aliases).__name__,
            "original_identifier_field_type":type(identifiers).__name__,
            "exact_match_in_2023_18_taxon_mapped_subset":matching,
            "source_2023_mapped_species_if_exact_match":original_tip_map.get(alias) if matching else None,
        })
    rows.sort(key=lambda r:r["biosample_accession"] or "")
    if {r["biosample_accession"] for r in rows}!=set(source_ids):
        raise ValueError("Original NCBI source aliases cannot be linked to verified BioSamples")
    return {
        "status":"NCBI_FIVE_SOURCE_ALIAS_DIAGNOSIS_ONLY",
        "n_fixed_source_ids":len(source_ids),
        "original_2023_mapped_tip_count":len(original_tip_map),
        "n_source_submitter_keys_available":sum(x["original_key_available"] for x in rows),
        "n_exact_old_2023_mapped_tips":sum(x["exact_match_in_2023_18_taxon_mapped_subset"] for x in rows),
        "source_alias_rows":rows,
        "no_sampledata_or_numeric_outcome_values_decoded":True,
        "not_full_2026_to_2023_genotype_identity":True,
    }


def request_json(url):
    req=urllib.request.Request(url,headers={
        "User-Agent":"TTF-Penstemon-original-sample-alias-source/1.0",
        "Accept":"application/json",
    })
    with urllib.request.urlopen(req,timeout=45) as resp:
        raw=resp.read(MAX_BYTES+1)
    if len(raw)>MAX_BYTES:raise ValueError("NCBI response exceeds source cap")
    return raw


def main():
    samples=source_samples()[:N]
    tipmap=original_2023_mapped_tips()
    term=" OR ".join(f"{acc}[Accession]" for acc in samples)
    q=urllib.parse.urlencode({"db":"biosample","term":term,"retmode":"json","retmax":"10"})
    try:source=request_json(NCBI_SEARCH+"?"+q)
    except urllib.error.HTTPError as err:
        if err.code in (400,401,403,404,429,503):
            print(json.dumps({"status":"NCBI_SOURCE_ALIAS_ESearch_UNAVAILABLE","http":err.code}))
            return
        raise
    uids=numeric_uids_from_esearch(source,samples)
    if uids["status"]!="NCBI_SOURCE_UIDS_RESOLVED":
        print(json.dumps({"status":"NCBI_SOURCE_ALIAS_ESearch_NOT_EVALUABLE","metadata":uids}))
        return
    q2=urllib.parse.urlencode({"db":"biosample","id":",".join(uids["ids"]),"retmode":"json"})
    try:raw=request_json(NCBI+"?"+q2)
    except urllib.error.HTTPError as err:
        if err.code in (400,401,403,404,429,503):
            print(json.dumps({"status":"NCBI_SOURCE_ALIAS_ESummary_UNAVAILABLE","http":err.code}))
            return
        raise
    joined=inspect_alias_fields(raw,samples,tipmap)
    print(json.dumps({"status":"PENSTEMON_SOURCE_ALIAS_PILOT_WITH_2023_ORIGINAL_METADATA",
        "NCBI_study":"PRJNA1479764","original_2023_groupmap_blob":SOURCE_BLOB,
        "source":joined,
        "all_52_2026_taxa_already_confirmed_on_separate_source_run":True,
        "no_flower_or_genotype_outcomes_read":True},indent=2,sort_keys=True))


if __name__=="__main__":main()
