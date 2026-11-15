"""Costus 2025 OSF scent × 2022 source-phylogeny source-key eligibility.

Original public OSF file downloads must match previously frozen OSF SHA256
before parsing. Reads only species-name columns and CSV headers/row counts;
never reads scent intensity, compound response, transcript expression or any
floral fitness response. Published chronogram is checksum-only and matches
Kay & Grossenbacher's 2022 Dryad chronogram byte-for-byte.
"""
from __future__ import annotations
import collections
import csv
import hashlib
import io
import json
import re
import urllib.error
import urllib.request

FILES={
 "2025_scent_samples":{
  "url":"https://osf.io/download/675c52794bafa93c71ea0a08/",
  "sha256":"161fbcee9d21dd2f4d0b9b62fe6e887f94c8bb190d22d17dd326fab1098d51b7",
  "size":7269,"kind":"species_csv",
 },
 "2025_floral_samples":{
  "url":"https://osf.io/download/675c5276cb373987a4de3063/",
  "sha256":"59097de1ba4328603fe0ee132df0f76033eaa94444ffd5485903ca9fb0db5c36",
  "size":5700,"kind":"species_csv",
 },
 "2025_species_mean_scent":{
  "url":"https://osf.io/download/675c526b379b48e27487bc42/",
  "sha256":"056b91173cd14907347557c267c802fa644c710830d039370b1b416e755f1af1",
  "size":6868,"kind":"species_csv",
 },
 "2025_scent_phylo_info":{
  "url":"https://osf.io/download/675c5278b0385fd0ba53bb5d/",
  "sha256":"9ce8f35b736652d3f0bbd2970707ce1aebebc0c4786962a23a983da7e81a8226",
  "size":1851,"kind":"species_csv",
 },
 "2025_gene_phylo_info":{
  "url":"https://osf.io/download/674f703d4afab2abfb9251a1/",
  "sha256":"375d2879f7c3217acacace997e2ecde9d084193c5251faf3d1eab8625f8f3243",
  "size":4687,"kind":"species_csv",
 },
 "2022_ancestral_chronogram_copy_in_2025_OSF":{
  "url":"https://osf.io/download/g3fve/",
  "sha256":"24f1cbfaa82ce4d33aad1ba80787cc611a6730aafb38f8244954c4d0d50ae7a1",
  "size":90309,"kind":"tree_checksum_only",
 },
}
SOURCE_2022_DRYAD_TREE_SHA="24f1cbfaa82ce4d33aad1ba80787cc611a6730aafb38f8244954c4d0d50ae7a1"
SPECIES_COLS={"species","species name","species_name","speciesname","species.id",
              "speciesid","taxon","taxa","sp","species.name","species_nm"}
INDEPENDENT_LABEL_COLS={"pollinator","pollinator_group","pollination","pollination.syndrome",
                        "pollination type","pollination_type","syndrome"}

def fetch_file(info):
    request=urllib.request.Request(info["url"],headers={
       "Accept":"application/octet-stream",
       "User-Agent":"TTF-Costus-OSF-original-plant-taxon-key-audit/1.0",
    })
    try:
        with urllib.request.urlopen(request,timeout=40) as response:
            data=response.read(info["size"]+1)
    except urllib.error.HTTPError as e:
        return None,{"status":"OSF_SOURCE_NOT_ACCESSIBLE","http":e.code,"url":info["url"]}
    if len(data)!=info["size"] or hashlib.sha256(data).hexdigest()!=info["sha256"]:
        raise ValueError("original OSF file size or SHA256 mismatch")
    return data,None

def normalize_taxon(s):
    return re.sub(r"\s+"," ",s.strip().replace("_"," ")).lower()

def inspect_csv_species_keys(raw):
    reader=csv.DictReader(io.StringIO(raw.decode("utf-8-sig")))
    hd=reader.fieldnames or []
    if not hd or len(hd)!=len(set(hd)):
        raise ValueError("duplicate or missing CSV header")
    fields={h.lower().strip():h for h in hd}
    species_field=next((fields[x] for x in SPECIES_COLS if x in fields),None)
    # deterministic decision on multiple aliases, never inferred by response.
    if not species_field:
        cand=[h for h in hd if re.search(r"species|^sp\.?$|^taxon$",h,re.I)]
        species_field=cand[0] if len(cand)==1 else None
    pollinator_headers=[h for h in hd if h.lower().strip() in INDEPENDENT_LABEL_COLS or "pollinator" in h.lower()]
    n=0
    unique=set()
    missing=0
    for row in reader:
        n+=1
        if species_field:
            val=(row.get(species_field) or "").strip()
            if val:unique.add(normalize_taxon(val))
            else:missing+=1
        # Explicitly DO NOT read other columns/values.
    return {
        "status":"SOURCE_SPECIES_KEYS_ONLY" if species_field else "SOURCE_HEADERS_ONLY_NO_UNAMBIGUOUS_SPECIES_COLUMN",
        "headers":hd,
        "row_count":n,
        "species_column":species_field,
        "unique_species_labels":sorted(unique) if species_field else None,
        "distinct_species_labels":len(unique) if species_field else None,
        "missing_species_values":missing if species_field else None,
        "pollinator_label_columns_PRESENT_NOT_INDEPENDENTLY_VERIFIED":pollinator_headers,
        "no_scent_or_pollination_values_read":True,
    }

def main():
    report={}
    for label,info in FILES.items():
        raw,problem=fetch_file(info)
        if problem:
            report[label]=problem
            continue
        out=inspect_csv_species_keys(raw) if info["kind"]=="species_csv" else {
            "status":"SOURCE_CHRONOGRAM_BYTES_VERIFIED_NO_TREE_RESPONSE_OR_TIP_VALUES_READ",
            "byte_identical_to_2022_original_dryad":info["sha256"]==SOURCE_2022_DRYAD_TREE_SHA,
            "no_scent_or_pollination_values_read":True,
        }
        report[label]={**out,"sha256":info["sha256"],"source_url":info["url"],"original_size":len(raw)}
    keyed={k:set(v.get("unique_species_labels") or []) for k,v in report.items()}
    candidates={k for k,v in keyed.items() if v}
    overlaps={
       "scent_samples_vs_means":sorted(keyed["2025_scent_samples"]&keyed["2025_species_mean_scent"]),
       "scent_samples_vs_scent_phylo_info":sorted(keyed["2025_scent_samples"]&keyed["2025_scent_phylo_info"]),
       "means_vs_scent_phylo_info":sorted(keyed["2025_species_mean_scent"]&keyed["2025_scent_phylo_info"]),
    }
    print(json.dumps({
       "status":"COSTUS_OSF_ORIGINAL_SPECIES_KEY_JOIN_ONLY",
       "sources":report,
       "species_keyed_source_count":len(candidates),
       "key_intersections":overlaps,
       "numbers_and_scent_responses_decoded":False,
       "pollinator_observation_provenance_not_verified":True,
       "phylogeny_topology_not_analyzed":True,
    },indent=2,sort_keys=True,ensure_ascii=False))

if __name__=="__main__":main()
