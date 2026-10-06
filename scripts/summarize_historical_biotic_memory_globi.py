#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import re
from collections import Counter, defaultdict
from pathlib import Path
from typing import Iterable


EXPECTED_RULE="ttf_historical_biotic_memory_globi_host_parser_rule_v0.1"
BAD_SECOND={"sp.","sp","spp.","spp","cf.","aff."}


def sha256_path(path: Path) -> str:
    h=hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda:handle.read(1<<20),b""):
            h.update(chunk)
    return h.hexdigest()


def read_species(path: Path) -> list[str]:
    rows=[x.strip() for x in path.read_text().splitlines() if x.strip()]
    if len(rows)!=len(set(rows)):
        raise RuntimeError("duplicate focal species")
    return rows


def read_host_families(path: Path) -> dict[str,str]:
    with path.open(newline="",encoding="utf-8") as handle:
        rows=list(csv.DictReader(handle))
    if set(rows[0])!={"species","host_family"}:
        raise RuntimeError("host-family schema drift")
    out={}
    for row in rows:
        sp=row["species"].strip()
        fam=row["host_family"].strip()
        if not sp or not fam or sp in out:
            raise RuntimeError("invalid host-family row")
        out[sp]=fam
    return out


def path_tokens(value: str) -> set[str]:
    return {x.strip() for x in str(value or "").split("|") if x.strip()}


def species_level_name(name: str) -> bool:
    parts=str(name or "").strip().split()
    if len(parts)<2:
        return False
    if parts[1].lower() in BAD_SECOND:
        return False
    if parts[0].lower().endswith("aceae"):
        return False
    return True


def iter_csv_files(root: Path) -> Iterable[Path]:
    yield from sorted(root.rglob("*.csv"))


def focal_from_filename(path: Path, species: list[str]) -> tuple[str,str]:
    m=re.match(r"^(\d{3})_(sourceTaxon|targetTaxon)\.csv$",path.name)
    if not m:
        raise RuntimeError(f"unexpected raw capture file name: {path.name}")
    idx=int(m.group(1))
    if not 0<=idx<len(species):
        raise RuntimeError("species index outside frozen panel")
    return species[idx],m.group(2)


def classify_row(
    row: dict[str,str],
    *,
    focal: str,
    family: str,
    rule: dict,
) -> tuple[str,str|None]:
    it=str(row.get("interaction_type","")).strip()
    source=str(row.get("source_taxon_name","")).strip()
    target=str(row.get("target_taxon_name","")).strip()

    forward=set(rule["accepted_orientations"]["focal_as_source"]["interaction_types"])
    inverse=set(rule["accepted_orientations"]["focal_as_target"]["interaction_types"])

    if source==focal and it in forward:
        plant=target
        plant_path=str(row.get("target_taxon_path",""))
    elif target==focal and it in inverse:
        plant=source
        plant_path=str(row.get("source_taxon_path",""))
    else:
        return "orientation_or_verb",None

    if family not in path_tokens(plant_path):
        return "host_family_mismatch",None
    if not species_level_name(plant):
        return "not_species_level_plant",None
    return "accepted",plant


def main() -> int:
    ap=argparse.ArgumentParser()
    ap.add_argument("--rule",type=Path,required=True)
    ap.add_argument("--species",type=Path,required=True)
    ap.add_argument("--host-families",type=Path,required=True)
    ap.add_argument("--raw-root",type=Path,required=True)
    ap.add_argument("--output",type=Path,required=True)
    args=ap.parse_args()

    rule=json.loads(args.rule.read_text())
    if rule.get("schema")!=EXPECTED_RULE:
        raise RuntimeError("unexpected parser rule")
    if not str(rule.get("status","")).startswith("FROZEN_"):
        raise RuntimeError("parser rule not frozen")
    if any(bool(v) for v in rule["response_firewall"].values()):
        raise RuntimeError("genetic response firewall is open")

    fi=rule["frozen_inputs"]
    if sha256_path(args.species)!=fi["species_sha256"]:
        raise RuntimeError("species SHA drift")
    if sha256_path(args.host_families)!=fi["host_family_sha256"]:
        raise RuntimeError("host-family SHA drift")
    species=read_species(args.species)
    families=read_host_families(args.host_families)
    if len(species)!=int(fi["species_count"]) or len(families)!=int(fi["host_family_rows"]):
        raise RuntimeError("frozen panel size drift")
    if set(species)!=set(families):
        raise RuntimeError("species/host-family membership mismatch")

    accepted=defaultdict(lambda:defaultdict(lambda:{"studies":set(),"interaction_types":set()}))
    rejected=defaultdict(Counter)
    raw_rows=Counter()
    files_seen=Counter()

    for path in iter_csv_files(args.raw_root):
        focal,direction=focal_from_filename(path,species)
        family=families[focal]
        files_seen[focal]+=1
        with path.open(newline="",encoding="utf-8-sig") as handle:
            reader=csv.DictReader(handle)
            required={
                "interaction_type","source_taxon_name","target_taxon_name",
                "source_taxon_path","target_taxon_path","study_title"
            }
            if not required.issubset(set(reader.fieldnames or [])):
                raise RuntimeError(f"GloBI schema drift in {path}")
            for row in reader:
                raw_rows[focal]+=1
                reason,plant=classify_row(row,focal=focal,family=family,rule=rule)
                if reason!="accepted":
                    rejected[focal][reason]+=1
                    continue
                rec=accepted[focal][str(plant)]
                rec["studies"].add(str(row.get("study_title","")).strip())
                rec["interaction_types"].add(str(row.get("interaction_type","")).strip())

    rows=[]
    for sp in species:
        hosts=[]
        for plant in sorted(accepted[sp]):
            rec=accepted[sp][plant]
            hosts.append({
                "raw_plant_name":plant,
                "studies":sorted(x for x in rec["studies"] if x),
                "interaction_types":sorted(rec["interaction_types"]),
            })
        rows.append({
            "species":sp,
            "host_family":families[sp],
            "raw_capture_files_seen":int(files_seen[sp]),
            "raw_rows":int(raw_rows[sp]),
            "accepted_species_level_host_count":len(hosts),
            "accepted_hosts":hosts,
            "rejected_reason_counts":dict(sorted(rejected[sp].items())),
        })

    captured=[x for x in rows if x["raw_capture_files_seen"]>0]
    payload={
        "schema":"ttf_historical_biotic_memory_globi_host_candidate_census_v0.1",
        "status":"EXTERNAL_HOST_CANDIDATE_CENSUS_COMPLETE",
        "captured_focal_species":len(captured),
        "full_panel_species":len(species),
        "focal_species_with_species_level_host":sum(x["accepted_species_level_host_count"]>0 for x in rows),
        "unique_focal_host_pairs":sum(x["accepted_species_level_host_count"] for x in rows),
        "rows":rows,
        "response_firewall":rule["response_firewall"],
    }
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(payload,indent=2,sort_keys=True)+"\n")
    print(json.dumps({
        "status":payload["status"],
        "captured_focal_species":payload["captured_focal_species"],
        "focal_species_with_species_level_host":payload["focal_species_with_species_level_host"],
        "unique_focal_host_pairs":payload["unique_focal_host_pairs"],
        "output_sha256":sha256_path(args.output),
    },sort_keys=True))
    return 0


if __name__=="__main__":
    raise SystemExit(main())
