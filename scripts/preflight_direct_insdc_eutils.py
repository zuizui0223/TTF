#!/usr/bin/env python3
"""Response-blind NCBI E-utilities preflight for the direct INSDC arm.

This step opens only the total number of records matching the frozen metadata query.
It does not fetch GenBank records, sequence identity, species counts, or geometry.
"""
from __future__ import annotations
import argparse, datetime as dt, json, urllib.parse, urllib.request

DEFAULT_QUERY = '(Animalia[Organism]) AND (COX1[Gene Name] OR CO1[Gene Name] OR COI[Gene Name] OR "cytochrome c oxidase subunit I"[Title])'
BASE = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi'


def esearch_count(term: str) -> int:
    data = urllib.parse.urlencode({'db':'nuccore','term':term,'retmode':'json','retmax':0}).encode()
    req = urllib.request.Request(BASE, data=data, headers={'User-Agent':'TTF-response-blind-source-preflight/0.1'})
    with urllib.request.urlopen(req, timeout=120) as r:
        payload = json.load(r)
    return int(payload['esearchresult']['count'])


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--query', default=DEFAULT_QUERY)
    ap.add_argument('--output', required=True)
    args=ap.parse_args()
    count=esearch_count(args.query)
    payload={
      'schema':'ttf_direct_insdc_eutils_query_preflight_v0.1',
      'retrieved_utc':dt.datetime.now(dt.timezone.utc).isoformat(),
      'query':args.query,
      'matching_record_count':count,
      'nucleotide_identity_opened':False,
      'species_geometry_opened':False,
      'empirical_ttf_opened':False,
      'next_step':'design_exhaustive_batched_metadata_fetch_from_exact_count'
    }
    with open(args.output,'w',encoding='utf-8') as f:
        json.dump(payload,f,indent=2,sort_keys=True); f.write('\n')

if __name__=='__main__': main()
