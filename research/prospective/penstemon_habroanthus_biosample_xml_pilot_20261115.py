"""Penstemon 2026 original BioSample XML identity availability PILOT.

Source BioProject PRJNA1479764: sample accessions obtained from the
verified ENA 52-run report. Select lexicographically FIRST FIVE distinct
BioSample IDs *before* reading any sample record. Query official ENA
Browser public XML only. No FASTQ, genotype, phenotype or trait data read.
"""
from __future__ import annotations
import collections
import csv
import io
import json
import urllib.error
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET

STUDY="PRJNA1479764"
LIMIT=5
RUN_URL="https://www.ebi.ac.uk/ena/portal/api/filereport"
SAMPLE_BASE="https://www.ebi.ac.uk/ena/browser/api/xml/"
MAX_XML=2_000_000


def get_source_samples()->list[str]:
    args=urllib.parse.urlencode({
      "accession":STUDY, "result":"read_run",
      "fields":"run_accession,sample_accession,study_accession",
      "format":"tsv",
    })
    req=urllib.request.Request(RUN_URL+"?"+args,headers={
      "User-Agent":"TTF-Penstemon-BioSample-identities/1.0"
    })
    with urllib.request.urlopen(req,timeout=40) as r:raw=r.read(250_001)
    if len(raw)>250_000:raise ValueError("ENA original run sample metadata too large")
    reader=csv.DictReader(io.StringIO(raw.decode("utf-8-sig")),delimiter="\t")
    if not {"run_accession","sample_accession","study_accession"}<=set(reader.fieldnames or []):
        raise ValueError("source run report headers differ")
    rows=list(reader)
    if len(rows)!=52 or len({r["run_accession"] for r in rows})!=52:
        raise ValueError("original ENA run accession coverage drift")
    if {r["study_accession"] for r in rows}!={STUDY}:
        raise ValueError("SRA BioProject differs")
    samples=sorted({r["sample_accession"] for r in rows})
    if len(samples)!=52 or any(not s for s in samples):
        raise ValueError("BioSample accession coverage differs")
    return samples


def inspect_xml(raw:bytes)->dict:
    if len(raw)>MAX_XML:raise ValueError("BioSample XML over safety cap")
    root=ET.fromstring(raw)
    # ENA Browser API XML sample records have SAMPLE or SAMPLE_SET roots.
    records=[root] if root.tag=="SAMPLE" else list(root.findall(".//SAMPLE"))
    if len(records)!=1:
        raise ValueError(f"unexpected BioSample XML record count: {len(records)}")
    item=records[0]
    sample_name=item.find("SAMPLE_NAME")
    names={}
    if sample_name is not None:
        for key in ("SCIENTIFIC_NAME","TAXON_ID","COMMON_NAME"):
            v=sample_name.findtext(key)
            if v:names[key]=v.strip()
    title=item.findtext("TITLE")
    # Only public voucher or submitter alias keys, not unrelated sample attributes.
    alias=item.attrib.get("alias")
    attrs=item.find("SAMPLE_ATTRIBUTES")
    attrnames=[]
    if attrs is not None:
        for row in attrs.findall("SAMPLE_ATTRIBUTE"):
            tag=row.findtext("TAG")
            if tag:attrnames.append(tag.strip())
    return {"source_xml_record_count":1,"scientific_name":names.get("SCIENTIFIC_NAME"),
      "taxon_id":names.get("TAXON_ID"),"title_present":bool(title),
      "submitter_alias_present":bool(alias),
      "available_attribute_field_names":sorted(set(attrnames)),
      "no_sequence_or_trait_values_read":True}


def main():
    allsamples=get_source_samples()
    chosen=allsamples[:LIMIT]
    records=[]
    for sample in chosen:
        url=SAMPLE_BASE+urllib.parse.quote(sample,safe="")
        req=urllib.request.Request(url,headers={"User-Agent":"TTF-Penstemon-ENA-BioSample-identity/1.0"})
        try:
            with urllib.request.urlopen(req,timeout=30) as response:
                data=response.read(MAX_XML+1)
        except urllib.error.HTTPError as err:
            records.append({"source_sample_accession":sample,
              "status":"PUBLIC_SAMPLE_XML_UNAVAILABLE",
              "http_status":err.code})
            continue
        try:
            body=inspect_xml(data)
        except (ET.ParseError,ValueError) as e:
            records.append({"source_sample_accession":sample,
                "status":"SOURCE_XML_SCHEMA_NOT_EVALUABLE",
                "error_type":type(e).__name__})
            continue
        records.append({"source_sample_accession":sample,"status":"PUBLIC_ORIGINAL_XML_IDENTITY",
                        **body})
    print(json.dumps({
       "status":"PENSTEMON_2026_BIOSAMPLE_XML_PILOT_SOURCE_ONLY",
       "original_source_bioproject":STUDY,"verified_source_samples":len(allsamples),
       "source_first_lexical_samples_attempted":len(chosen),
       "sample_xml_source_results":records,
       "no_biology_outcomes_read":True,
       "no_genotype_or_reads_downloaded":True,
       "2023_older_tree_specimen_join_not_yet_authenticated":True,
    },sort_keys=True,indent=2))

if __name__=="__main__":main()
