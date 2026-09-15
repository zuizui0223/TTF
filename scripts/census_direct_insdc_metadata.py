#!/usr/bin/env python3
"""Build a metadata-only census from GenBank flatfiles.

No ORIGIN sequence content is parsed or persisted. This script is source-feasibility
infrastructure only; it does not authorize genetic-distance or TTF inference.
"""
from __future__ import annotations
import argparse, gzip, json, math, re
from collections import defaultdict
from pathlib import Path
from ttf.direct_insdc_metadata import iter_genbank_metadata

LATLON = re.compile(r'^\s*([+-]?\d+(?:\.\d+)?)\s*([NS])?[, ]+\s*([+-]?\d+(?:\.\d+)?)\s*([EW])?\s*$')


def parse_latlon(s: str):
    m = LATLON.match(s)
    if not m:
        return None
    lat, ns, lon, ew = m.groups(); lat=float(lat); lon=float(lon)
    if ns == 'S': lat=-abs(lat)
    if ns == 'N': lat=abs(lat)
    if ew == 'W': lon=-abs(lon)
    if ew == 'E': lon=abs(lon)
    if not (-90 <= lat <= 90 and -180 <= lon <= 180): return None
    return (round(lat, 6), round(lon, 6))


def open_text(path: Path):
    return gzip.open(path, 'rt', encoding='utf-8', errors='replace') if path.suffix == '.gz' else path.open('rt', encoding='utf-8', errors='replace')


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('flatfiles', nargs='+', type=Path)
    ap.add_argument('--exclude-species', type=Path)
    ap.add_argument('--min-localities', type=int, default=12)
    ap.add_argument('--output', type=Path, required=True)
    args=ap.parse_args()
    excluded=set()
    if args.exclude_species:
        excluded={x.strip() for x in args.exclude_species.read_text().splitlines() if x.strip()}
    locs=defaultdict(set); records=0; coi=0
    for path in args.flatfiles:
        with open_text(path) as fh:
            for rec in iter_genbank_metadata(fh):
                records += 1
                if not rec.coi_family_annotated: continue
                coi += 1
                if 'Animalia' not in rec.taxonomy or 'Aves' in rec.taxonomy or 'Chiroptera' in rec.taxonomy: continue
                if rec.organism in excluded: continue
                xy=parse_latlon(rec.lat_lon)
                if xy is not None: locs[rec.organism].add(xy)
    counts=sorted((sp,len(v)) for sp,v in locs.items())
    eligible=[(sp,n) for sp,n in counts if n >= args.min_localities]
    payload={
      'schema':'ttf_direct_insdc_metadata_census_v0.1',
      'sequence_identity_opened':False,
      'empirical_genetic_outcome_opened':False,
      'records_with_required_metadata':records,
      'coi_family_records':coi,
      'species_with_georeferenced_coi':len(counts),
      'species_meeting_locality_floor':len(eligible),
      'minimum_localities':args.min_localities,
      'decision':'PASS_TO_GEOMETRY_QUALIFICATION' if len(eligible) >= 150 else 'NOT_EVALUABLE_SOURCE_FEASIBILITY',
      'species_locality_counts':eligible,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2, sort_keys=True)+'\n')

if __name__ == '__main__': main()
