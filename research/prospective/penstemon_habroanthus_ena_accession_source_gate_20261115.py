"""Penstemon Habroanthus 2026 genome BioProject source-only accession gate.

Source: Stone et al. 2026 bioRxiv DOI 10.64898/2026.07.21.739878,
stating NCBI BioProject PRJNA1479764, 52 sampled accessions (20
Habroanthus species + one outgroup). Query public ENA run metadata ONLY,
never downloading reads, genotypes or genetic outcomes.
"""
from __future__ import annotations
import collections
import csv
import io
import json
import urllib.error
import urllib.parse
import urllib.request

BIOPROJECT="PRJNA1479764"
ENA="https://www.ebi.ac.uk/ena/portal/api/filereport"
FIELDS=("run_accession","sample_accession","sample_alias","scientific_name",
        "study_accession","experiment_accession")
MANUSCRIPT="10.64898/2026.07.21.739878"
MAX_REPORT_BYTES=2_000_000


def parse_public_metadata(raw:bytes)->dict:
    text=raw.decode("utf-8-sig")
    rows=list(csv.DictReader(io.StringIO(text),delimiter="\t"))
    if not rows:
        return {"status":"NO_PUBLIC_RUN_METADATA","run_count":0,
                "sample_accessions":0,"source_species_labels":[],
                "not_genome_absence":True}
    seen=set()
    counts=collections.Counter()
    species=collections.Counter()
    sample_aliases=collections.Counter()
    studyids=set()
    for row in rows:
        run=(row.get("run_accession") or "").strip()
        sample=(row.get("sample_accession") or "").strip()
        study=(row.get("study_accession") or "").strip()
        name=(row.get("scientific_name") or "").strip()
        alias=(row.get("sample_alias") or "").strip()
        if not run or not sample or run in seen:
            raise ValueError("Duplicate or missing original ENA run/sample accession")
        seen.add(run)
        counts[sample]+=1
        if study:studyids.add(study)
        if name:species[name]+=1
        if alias:sample_aliases[alias]+=1
    return {
        "status":"ENA_PUBLIC_BIOPROJECT_SOURCE_METADATA_ONLY",
        "run_count":len(rows),
        "sample_accession_count":len(counts),
        "source_species_label_count":len(species),
        "source_species_labels":dict(sorted(species.items())),
        "study_accessions":sorted(studyids),
        "sample_alias_count":len(sample_aliases),
        "sample_alias_duplicates":sorted([k for k,v in sample_aliases.items() if v>1]),
        "source_field_names":list(rows[0]),
        "run_to_sample_original_crosswalk_not_exported":True,
        "published_52_accessions_not_assumed_to_equal_run_count":True,
        "no_biological_outcomes_read":True,
    }


def main():
    params=urllib.parse.urlencode({
      "accession":BIOPROJECT,"result":"read_run",
      "fields":",".join(FIELDS),"format":"tsv"
    })
    url=ENA+"?"+params
    req=urllib.request.Request(url,headers={"Accept":"text/tab-separated-values",
               "User-Agent":"TTF-Penstemon-source-metadata/1.0"})
    try:
        with urllib.request.urlopen(req,timeout=40) as resp:
            data=resp.read(MAX_REPORT_BYTES+1)
    except urllib.error.HTTPError as err:
        if err.code in (400,401,403,404,429):
            print(json.dumps({
                "status":"ENA_PUBLIC_METADATA_REQUEST_BLOCKED",
                "http_status":err.code,
                "BioProject":BIOPROJECT,
                "fields":FIELDS,
                "does_not_establish_52_missing_accessions":True,
                "no_biological_outcomes_read":True,
            },indent=2))
            return
        raise
    if len(data)>MAX_REPORT_BYTES:raise ValueError("ENA response unexpectedly large")
    try:
        parsed=parse_public_metadata(data)
    except UnicodeDecodeError:
        raise ValueError("not ENA UTF-8 source metadata")
    print(json.dumps({
      "source_bioproject":BIOPROJECT,
      "source_manuscript_doi":MANUSCRIPT,
      "reported_author_accessions":52,
      "reported_Habroanthus_species":20,
      "reported_outgroup_species":"Penstemon palmeri",
      "source_ena_query_fields":FIELDS,
      "result":parsed,
      "no_sequences_downloaded":True,
    },indent=2,sort_keys=True,ensure_ascii=False))


if __name__=="__main__":main()
