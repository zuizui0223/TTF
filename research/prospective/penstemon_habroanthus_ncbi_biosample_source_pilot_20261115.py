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
    q=urllib.parse.urlencode({
       "db":"biosample","id":",".join(samples),"retmode":"json",
       "tool":"TTF_CREST_habroanthus_source_pilot",
       "email":"", # Anonymous public E-utility rate limit, no API key.
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
         "result":result,
         "no_sequence_genotype_or_flower_outcomes_opened":True,
         "not_2023_specimen_join":True,
    },indent=2,sort_keys=True))


if __name__=="__main__":main()
