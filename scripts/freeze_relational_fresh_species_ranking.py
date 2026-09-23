#!/usr/bin/env python3
"""Build a response-blind fresh-species ranking for Relational TTF.

This script intentionally operates only on species names and exclusion lists.
It never reads sequence identity, genetic distance, or T_st.
"""
from __future__ import annotations
import argparse, csv, hashlib, json
from pathlib import Path

TAG = "relational-ttf-v0.1"

def rank_key(source_sha: str, species: str) -> tuple[str, str]:
    token = f"{TAG}|{source_sha}|{species}".encode()
    return hashlib.sha256(token).hexdigest(), species

def read_species(path: Path) -> set[str]:
    if not path.exists():
        return set()
    if path.suffix.lower() == ".json":
        obj=json.loads(path.read_text())
        if isinstance(obj, list): return {str(x).strip() for x in obj if str(x).strip()}
        for key in ("species","selected_species","excluded_species"):
            if isinstance(obj.get(key), list):
                return {str(x).strip() for x in obj[key] if str(x).strip()}
        raise ValueError(f"no species list in {path}")
    with path.open(newline="") as f:
        rows=list(csv.DictReader(f))
    if not rows or "species" not in rows[0]:
        raise ValueError(f"{path} must contain a species column")
    return {r["species"].strip() for r in rows if r["species"].strip()}

def main() -> None:
    p=argparse.ArgumentParser()
    p.add_argument("--universe", type=Path, required=True, help="response-blind CSV/JSON species universe")
    p.add_argument("--source-sha256", required=True)
    p.add_argument("--exclude", type=Path, action="append", default=[])
    p.add_argument("--n", type=int, default=1000)
    p.add_argument("--output", type=Path, required=True)
    a=p.parse_args()
    universe=read_species(a.universe)
    excluded=set()
    for path in a.exclude: excluded |= read_species(path)
    eligible=sorted(universe-excluded, key=lambda x: rank_key(a.source_sha256,x))
    chosen=eligible[:a.n]
    payload={
      "schema":"ttf_relational_fresh_species_ranking_v0.1",
      "tag":TAG,
      "source_sha256":a.source_sha256,
      "universe_n":len(universe),
      "excluded_n":len(universe & excluded),
      "eligible_n":len(eligible),
      "requested_n":a.n,
      "selected_n":len(chosen),
      "ranking_rule":"SHA256('relational-ttf-v0.1|<source_sha256>|<species>'), species tie-break",
      "selected_species":chosen,
      "selected_species_digest_sha256":hashlib.sha256(("\n".join(chosen)+"\n").encode()).hexdigest(),
      "firewall":{"genetic_response_opened":False,"T_st_computed":False}
    }
    a.output.parent.mkdir(parents=True,exist_ok=True)
    a.output.write_text(json.dumps(payload,indent=2)+"\n")

if __name__=="__main__":
    main()
