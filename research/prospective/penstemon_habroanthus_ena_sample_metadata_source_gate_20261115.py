"""ENA sample-level taxon/voucher metadata for Penstemon Habroanthus BioProject.

Source BioProject PRJNA1479764 published in Stone et al 2026. Query only
documented ENA public portal fields; no sequences/genotypes are requested.
The metadata source might lack scientific names/voucher identifiers, in which
case the join to a published phenotype accession is NOT_EVALUABLE.
"""
from __future__ import annotations
import csv
import io
import json
import urllib.error
import urllib.parse
import urllib.request
from collections import Counter

BIOPROJECT="PRJNA1479764"
BASE="https://www.ebi.ac.uk/ena/portal/api"
DESIRED=("sample_accession","sample_alias","sample_title","scientific_name",
         "tax_id","sample_description","study_accession","collection_date")
ALWAYS=("sample_accession","scientific_name")
MAX_BYTES=3_000_000

def fetch_json(url):
    req=urllib.request.Request(url,headers={"User-Agent":"TTF-Penstemon-sample-source-gate/1.0"})
    with urllib.request.urlopen(req,timeout=45) as resp:raw=resp.read(MAX_BYTES+1)
    if len(raw)>MAX_BYTES:raise ValueError("ENA field dictionary too large")
    value=raw.decode("utf-8-sig").strip()
    try:
        return json.loads(value)
    except json.JSONDecodeError:
        # ENA returnFields is also served as a tabular plain-text field list.
        # Do not guess biological records from malformed response bodies.
        if value.lower().startswith(("<html","<!doctype")):
            raise ValueError("ENA returned HTML instead of public source field dictionary")
        lines=value.splitlines()
        if not 2<=len(lines)<=3000:
            raise ValueError("Unrecognized ENA returnFields text response shape")
        first=lines[0].lower()
        if not any(token in first for token in ("name","field","description")):
            raise ValueError("ENA returnFields text has no recognizable header")
        rd=csv.DictReader(io.StringIO(value),delimiter="\t")
        names=[]
        for row in rd:
            candidate=next((v for k,v in row.items()
                           if k and k.lower().strip() in
                           ("name","field","fieldname","field name","columnid")),None)
            if candidate and candidate.strip():
                names.append({"name":candidate.strip()})
        if not names:
            raise ValueError("No ENA field names could be extracted from tabular response")
        return names

def discover_fields(payload):
    if isinstance(payload,list):
        names={x.get("name") or x.get("field") for x in payload if isinstance(x,dict)}
    elif isinstance(payload,dict):
        names=set(payload)
    else:
        raise ValueError("Unknown ENA sample field schema")
    if not set(ALWAYS)<=names:
        raise ValueError(f"Essential ENA sample identity fields missing: {sorted(set(ALWAYS)-names)}")
    return [f for f in DESIRED if f in names]

def count_key_fields(raw):
    if not raw or raw.startswith(b"<html") or raw.startswith(b"<!DOCTYPE"):
        raise ValueError("ENA sample metadata API returned HTML not TSV")
    reader=csv.DictReader(io.StringIO(raw.decode("utf-8-sig")),delimiter="\t")
    headers=reader.fieldnames or []
    if "sample_accession" not in headers or len(headers)!=len(set(headers)):
        raise ValueError("Sample metadata headers missing or duplicated")
    ids=set()
    present=Counter()
    taxa=Counter()
    for row in reader:
        acc=(row.get("sample_accession") or "").strip()
        if not acc or acc in ids:raise ValueError("Duplicate or absent ENA BioSample accession")
        ids.add(acc)
        for f in headers:
            val=(row.get(f) or "").strip()
            if val:present[f]+=1
        sci=(row.get("scientific_name") or "").strip()
        if sci:taxa[sci]+=1
    return {"n_unique_samples":len(ids),"source_fields":headers,
            "source_nonempty_counts":dict(sorted(present.items())),
            "original_scientific_name_labels":dict(sorted(taxa.items())),
            "sample_accession_values_not_exported":True,
            "original_voucher_join_not_verified":True,
            "not_a_floral_trait_measurement":True}

def main():
    url=f"{BASE}/returnFields?result=sample"
    try:
        fields=discover_fields(fetch_json(url))
    except (urllib.error.HTTPError,ValueError,KeyError) as err:
        print(json.dumps({"status":"ENA_SAMPLE_FIELD_SCHEMA_NOT_ACCESSIBLE",
         "error_type":type(err).__name__,
         "no_negative_biological_claim":True},indent=2))
        return
    query=urllib.parse.urlencode({
        "accession":BIOPROJECT,"result":"sample",
        "fields":",".join(fields),"format":"tsv"
    })
    target=f"{BASE}/filereport?{query}"
    try:
        request=urllib.request.Request(target,headers={"User-Agent":"TTF-Penstemon-sample-source-gate/1.0"})
        with urllib.request.urlopen(request,timeout=50) as r:raw=r.read(MAX_BYTES+1)
    except urllib.error.HTTPError as err:
        if err.code in (400,401,403,404,429):
            print(json.dumps({"status":"ENA_SAMPLE_FILE_REPORT_BLOCKED",
              "http":err.code,"fields":fields,
              "not_absence_of_bioproject":True},indent=2))
            return
        raise
    if len(raw)>MAX_BYTES:raise ValueError("ENA sample source too large")
    result=count_key_fields(raw)
    print(json.dumps({"status":"ENA_HABROANTHUS_SAMPLE_SOURCE_METADATA_ONLY",
         "source_bioproject":BIOPROJECT,
         "source_requested_fields":fields,"result":result,
         "no_sequence_or_genotype_or_floral_outcome_read":True},indent=2,sort_keys=True))

if __name__=="__main__":main()
