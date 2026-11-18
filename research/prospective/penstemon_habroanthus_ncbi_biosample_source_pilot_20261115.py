"""Penstemon 2026 original NCBI BioSample accession-identity PILOT.

Take the lexicographically first five *validated* PRJNA1479764 BioSample
accessions from the independent ENA read_run report. Query NCBI's officially
documented E-utilities biosample ESummary by accession ID; parse only taxon
name / accession and the presence of voucher/alias attributes (NOT values).
If absent or unsupported, stop as accession-identity NOT_EVALUABLE; no
phenotype/genotype data or sequence outcomes are opened.
"""
from __future__ import annotations
import csv
import io
import json
import urllib.error
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET

STUDY="PRJNA1479764"
LIMIT=5
SOURCE="https://www.ebi.ac.uk/ena/portal/api/filereport"
NCBI_SEARCH="https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
NCBI="https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"
FIELDS=("run_accession","sample_accession","study_accession")
MAX_BYTES=2_000_000


def source_samples():
    query=urllib.parse.urlencode({"accession":STUDY,"result":"read_run",
       "fields":",".join(FIELDS),"format":"tsv"})
    req=urllib.request.Request(SOURCE+"?"+query,headers={"User-Agent":"TTF-Habroanthus-ncbi-metadata/1.0"})
    with urllib.request.urlopen(req,timeout=40) as response:
        raw=response.read(250_001)
    if len(raw)>250_000:raise ValueError("ENA source run metadata too large")
    dr=csv.DictReader(io.StringIO(raw.decode("utf-8-sig")),delimiter="\t")
    if set(FIELDS)-set(dr.fieldnames or []):raise ValueError("Original ENA fields differ")
    rows=list(dr)
    if len(rows)!=52 or len({r["run_accession"] for r in rows})!=52:
        raise ValueError("Original 52 run source count drift")
    if {r["study_accession"] for r in rows}!={STUDY}:
        raise ValueError("BioProject authority mismatch")
    ids=sorted({r["sample_accession"] for r in rows})
    if len(ids)!=52 or any(not a.startswith("SAMN") for a in ids):
        raise ValueError("BioSample accession type/count mismatch")
    return ids


def numeric_uids_from_esearch(raw:bytes,source_accessions:list[str])->dict:
    """Resolve SAMN source accessions to numeric NCBI Entrez UIDs, before ESummary."""
    try: data=json.loads(raw)
    except json.JSONDecodeError:
        raise ValueError("ESearch source metadata is not JSON")
    res=data.get("esearchresult") or {}
    if not isinstance(res,dict):raise ValueError("ESearch result changed")
    if "errorlist" in res:
        return {"status":"NCBI_ESEARCH_SOURCE_FIELD_NOT_EVALUABLE",
                "error_keys":sorted(res["errorlist"])}
    ids=res.get("idlist") or []
    try: count=int(res.get("count","-1"))
    except (TypeError,ValueError):raise ValueError("Invalid ESearch count")
    if count!=len(source_accessions) or len(ids)!=len(source_accessions):
        return {"status":"NCBI_ESEARCH_UID_COVERAGE_NOT_EVALUABLE",
                "source_accessions":len(source_accessions),
                "esearch_count":count,"numeric_uids_found":len(ids),
                "not_a_biological_null":True}
    if len(set(ids))!=len(ids) or not all(str(x).isdigit() for x in ids):
        raise ValueError("ESearch UID result not unique or numeric")
    return {"status":"NCBI_SOURCE_UIDS_RESOLVED",
            "ids":[str(x) for x in ids]}


def source_metadata_aggregate(raw:bytes,attempted:list[str])->dict:
    if len(raw)>MAX_BYTES:raise ValueError("NCBI ESummary metadata exceeded cap")
    try:obj=json.loads(raw)
    except json.JSONDecodeError:
        raise ValueError("NCBI ESummary returned non-JSON metadata")
    if not isinstance(obj,dict):
        raise ValueError("NCBI ESummary top-level schema not object")
    if "error" in obj:
        return {"status":"NCBI_ESUMMARY_REPORTED_ERROR",
                "message_excerpt":str(obj["error"])[:180]}
    d=obj.get("result") or {}
    if not isinstance(d,dict):
        raise ValueError("NCBI ESummary result schema changed")
    uids=d.get("uids") or []
    if not isinstance(uids,list):
        raise ValueError("NCBI returned non-list BioSample identifiers")
    records=[]
    reported=set()
    for uid in uids:
        entry=d.get(str(uid)) or {}
        if not isinstance(entry,dict):
            raise ValueError("NCBI BioSample entry not object")
        accession=entry.get("accession") or entry.get("Accession")
        if accession:
            reported.add(accession)
        organism=entry.get("organism") or entry.get("Organism")
        if isinstance(organism,dict):
            taxon=organism.get("scientificname") or organism.get("taxonomyname") or organism.get("name")
        else:
            taxon=organism if isinstance(organism,str) else None
        records.append({"source_biosample_accession":accession,
                "scientific_name":taxon,"has_submission_title":bool(entry.get("title")),
                "metadata_keys":sorted(entry.keys()),
                "raw_sample_attribute_values_not_decoded":True})
    extras=sorted(reported-set(attempted))
    if extras:
        raise ValueError("NCBI ESummary returned accession outside fixed five source IDs")
    return {"status":"NCBI_BIOSAMPLE_SOURCE_SUMMARY_ONLY",
            "requested":len(attempted),"result_records":len(records),
            "accessions_with_taxon_identity":sum(bool(r["scientific_name"]) for r in records),
            "records":records,
            "unmatched_source_accessions":sorted(set(attempted)-reported),
            "no_sequence_or_genotype_or_phenotype_values_read":True}


def main():
    samples=source_samples()[:LIMIT]
    term=" OR ".join(f"{sample}[Accession]" for sample in samples)
    searchquery=urllib.parse.urlencode({
        "db":"biosample","term":term,"retmode":"json",
        "retmax":"10","tool":"TTF_CREST_habroanthus_source_pilot",
    })
    searchreq=urllib.request.Request(NCBI_SEARCH+"?"+searchquery,
      headers={"User-Agent":"TTF-Habroanthus-BioSample-ESearch/1.0","Accept":"application/json"})
    try:
        with urllib.request.urlopen(searchreq,timeout=45) as response:
            searchraw=response.read(MAX_BYTES+1)
    except urllib.error.HTTPError as e:
        if e.code in (400,401,403,404,429,503):
            print(json.dumps({"status":"NCBI_ESearch_SOURCE_ACCESS_BLOCKED",
                "http_status":e.code,"source_samples_confirmed":52,
                "not_a_biological_null":True},indent=2))
            return
        raise
    resolved=numeric_uids_from_esearch(searchraw,samples)
    if resolved["status"]!="NCBI_SOURCE_UIDS_RESOLVED":
        print(json.dumps({"status":"NCBI_BIOSAMPLE_UID_JOIN_NOT_EVALUABLE",
            "original_52_samples_verified":True,
            "esearch":resolved,"no_biological_outcomes_opened":True},indent=2))
        return
    q=urllib.parse.urlencode({
        "db":"biosample","id":",".join(resolved["ids"]),"retmode":"json",
        "tool":"TTF_CREST_habroanthus_source_pilot",
    })
    req=urllib.request.Request(NCBI+"?"+q,headers={
        "User-Agent":"TTF-Habroanthus-public-biosample-source/1.0",
        "Accept":"application/json",
    })
    try:
        with urllib.request.urlopen(req,timeout=45) as resp:
            raw=resp.read(MAX_BYTES+1)
    except urllib.error.HTTPError as err:
        if err.code in (400,401,403,404,429,503):
            print(json.dumps({"status":"NCBI_BIOSAMPLE_ORIGINAL_SUMMARY_ACCESS_BLOCKED",
                "http_status":err.code,
                "source_samples_confirmed":52,"first_fixed_sample_count":LIMIT,
                "not_a_missing_genome_biological_result":True},indent=2))
            return
        raise
    result=source_metadata_aggregate(raw,samples)
    print(json.dumps({
         "status":"PENSTEMON_NCBI_FIVE_SOURCE_BIOSAMPLE_METADATA_PILOT",
         "BioProject":STUDY,"verified_ena_source_samples":52,
         "original_source_first_five_samples":samples,
         "NCBI_source_numeric_UIDs_resolved":len(resolved["ids"]),
         "result":result,
         "no_sequence_genotype_or_flower_outcomes_opened":True,
         "not_2023_specimen_join":True,
    },indent=2,sort_keys=True))


if __name__=="__main__":main()
