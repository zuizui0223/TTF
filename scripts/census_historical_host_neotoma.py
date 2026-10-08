#!/usr/bin/env python3
"""Response-blind Neotoma LGM coverage census for frozen host plants.

This script reads only the frozen insect-host panel and public paleoecological
occurrence metadata. It never reads nucleotide identity, genetic distances,
post-IBD turnover, or any ecology-genetics association.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import os
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from typing import Any

import requests

API = "https://api.neotomadb.org/v2.0/data/occurrences"
PAGE_LIMIT = 10000


def checkpoint_name(genus: str, dataset_type: str) -> str:
    """Content-independent deterministic key for an exact frozen query."""
    key=json.dumps([str(genus),str(dataset_type)],separators=(",",":"))
    return hashlib.sha256(key.encode("utf-8")).hexdigest()+".json"


def read_completed_checkpoint(
    checkpoint_dir: Path, genus: str, dataset_type: str, age_young: int, age_old: int
) -> dict[str, Any] | None:
    path=Path(checkpoint_dir)/checkpoint_name(genus,dataset_type)
    if not path.exists():
        return None
    rec=json.loads(path.read_text(encoding="utf-8"))
    expected={
        "schema":"ttf_historical_host_neotoma_query_checkpoint_v0.1",
        "genus":genus,
        "dataset_type":dataset_type,
        "age_young":int(age_young),
        "age_old":int(age_old),
    }
    for key,val in expected.items():
        if rec.get(key)!=val:
            raise RuntimeError(f"Neotoma checkpoint source query drift: {path.name} {key}")
    data=rec.get("result")
    if not isinstance(data,dict) or data.get("genus")!=genus or data.get("dataset_type")!=dataset_type:
        raise RuntimeError("Neotoma checkpoint result identity drift")
    if not isinstance(data.get("pairs"),list):
        raise RuntimeError("Neotoma checkpoint pairs must be a list")
    return data


def write_completed_checkpoint(
    checkpoint_dir: Path, genus: str, dataset_type: str,
    age_young: int, age_old: int, result: dict[str,Any]
) -> Path:
    if result.get("genus")!=genus or result.get("dataset_type")!=dataset_type:
        raise RuntimeError("Neotoma checkpoint result does not match requested source query")
    checkpoint_dir=Path(checkpoint_dir)
    checkpoint_dir.mkdir(parents=True,exist_ok=True)
    path=checkpoint_dir/checkpoint_name(genus,dataset_type)
    payload={
        "schema":"ttf_historical_host_neotoma_query_checkpoint_v0.1",
        "genus":genus,
        "dataset_type":dataset_type,
        "age_young":int(age_young),
        "age_old":int(age_old),
        "result":result,
    }
    temp=path.with_suffix(".tmp")
    temp.write_text(json.dumps(payload,indent=2,sort_keys=True)+"\n",encoding="utf-8")
    os.replace(temp,path)
    return path


def sha256_path(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def norm_name(x: str) -> str:
    return " ".join(str(x or "").split())


def api_dataset_type(label: str) -> str:
    label = norm_name(label).lower()
    if label == "plant macrofossils":
        return "plant macrofossil"
    return label


def parse_records(payload: dict[str, Any]) -> list[dict[str, Any]]:
    data = payload.get("data", [])
    if not isinstance(data, list):
        raise RuntimeError("Neotoma response data is not a list")
    out: list[dict[str, Any]] = []
    for row in data:
        if not isinstance(row, dict):
            continue
        site = row.get("site") if isinstance(row.get("site"), dict) else {}
        taxon = row.get("taxon") if isinstance(row.get("taxon"), dict) else {}
        sample = row.get("sample") if isinstance(row.get("sample"), dict) else {}
        ages = row.get("ages") if isinstance(row.get("ages"), dict) else {}
        taxonname = norm_name(taxon.get("taxonname") or sample.get("taxonname") or "")
        siteid = site.get("siteid")
        if not taxonname or siteid is None:
            continue
        out.append(
            {
                "taxonname": taxonname,
                "siteid": int(siteid),
                "datasettype": norm_name(site.get("datasettype", "")),
                "age": ages.get("age"),
                "ageolder": ages.get("ageolder"),
                "ageyounger": ages.get("ageyounger"),
            }
        )
    return out


def get_json(session: requests.Session, params: dict[str, Any], *, retries: int = 4) -> dict[str, Any]:
    last: Exception | None = None
    for attempt in range(retries):
        try:
            resp = session.get(API, params=params, timeout=90)
            resp.raise_for_status()
            payload = resp.json()
            if not isinstance(payload, dict):
                raise RuntimeError("Neotoma response is not an object")
            return payload
        except Exception as exc:  # pragma: no cover - network path
            last = exc
            if attempt + 1 < retries:
                time.sleep(2 ** attempt)
    raise RuntimeError(f"Neotoma request failed after {retries} attempts: {last}")


def fetch_genus(
    genus: str,
    dataset_type: str,
    age_young: int,
    age_old: int,
) -> dict[str, Any]:
    session = requests.Session()
    session.headers.update({"User-Agent": "TTF-historical-host-connectivity/0.1"})
    offset = 0
    records: list[dict[str, Any]] = []
    pages = 0
    while True:
        params = {
            "taxonname": f"{genus}%",
            "datasettype": dataset_type,
            "ageyoung": int(age_young),
            "ageold": int(age_old),
            "limit": PAGE_LIMIT,
            "offset": offset,
        }
        payload = get_json(session, params)
        page = parse_records(payload)
        records.extend(page)
        pages += 1
        raw_data = payload.get("data", [])
        if not isinstance(raw_data, list) or len(raw_data) < PAGE_LIMIT:
            break
        offset += PAGE_LIMIT
        if pages >= 25:
            raise RuntimeError(f"Neotoma pagination ceiling reached for {genus} / {dataset_type}")
    unique = {
        (r["taxonname"], int(r["siteid"]))
        for r in records
        if norm_name(r["taxonname"]).lower().startswith(genus.lower())
    }
    return {
        "genus": genus,
        "dataset_type": dataset_type,
        "records": len(records),
        "pages": pages,
        "pairs": sorted(unique),
    }


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--panel", type=Path, required=True)
    ap.add_argument("--rule", type=Path, required=True)
    ap.add_argument("--identity-correction", type=Path, required=True)
    ap.add_argument("--wcvp-names", type=Path, required=True)
    ap.add_argument("--output-dir", type=Path, required=True)
    ap.add_argument("--workers", type=int, default=6)
    ap.add_argument("--checkpoint-dir", type=Path, default=None)
    args = ap.parse_args()

    rule = json.loads(args.rule.read_text())
    if rule.get("schema") not in {"ttf_genetic_historical_host_neotoma_validation_rule_v0.1", "ttf_genetic_historical_host_neotoma_validation_rule_v0.2"}:
        raise RuntimeError("unexpected Neotoma validation rule")
    if rule.get("status") not in {"FROZEN_BEFORE_ANY_NEOTOMA_HOST_COVERAGE_RESULT_OR_PALAEO_HOST_RESISTANCE_RESULT", "FROZEN_BEFORE_V0_2_NEOTOMA_EXTERNAL_QUERY"}:
        raise RuntimeError("Neotoma rule is not prospectively frozen")

    with args.panel.open(newline="", encoding="utf-8") as f:
        panel = list(csv.DictReader(f))
    correction = json.loads(args.identity_correction.read_text())
    if correction.get("schema") not in {"ttf_genetic_historical_host_panel_identity_correction_v0.1", "ttf_genetic_historical_host_panel_identity_binding_v0.2"}:
        raise RuntimeError("unexpected panel identity correction")
    expected = rule["primary_panel"]
    ident = correction["scientific_identity"]
    if len(panel) != int(expected["species"]) or len(panel) != int(ident["insects"]):
        raise RuntimeError("frozen panel species count drift")

    if rule.get("schema", "").endswith("v0.2") and sha256_path(args.panel) != expected["sha256"]:
        raise RuntimeError("v0.2 panel byte SHA drift")

    species = sorted(norm_name(row["species"]) for row in panel)
    species_digest = hashlib.sha256(("\n".join(species)+"\n").encode()).hexdigest()
    if species_digest != ident["species_digest_sha256"]:
        raise RuntimeError("frozen panel species identity drift")

    logical = []
    panel_ids: dict[str, list[str]] = {}
    for row in panel:
        insect = norm_name(row["species"])
        ids = sorted(x.strip() for x in row["accepted_host_ids"].split(";") if x.strip())
        if int(row["n_hosts"]) != len(ids):
            raise RuntimeError(f"host ID count drift for {insect}")
        panel_ids[insect] = ids
        logical.append(f"{insect}|{';'.join(ids)}")
    logical_digest = hashlib.sha256(("\n".join(sorted(logical))+"\n").encode()).hexdigest()
    if logical_digest != ident["logical_identity_digest_sha256"]:
        raise RuntimeError("frozen insect x host-ID identity drift")

    with args.wcvp_names.open(newline="", encoding="utf-8") as f:
        name_rows = list(csv.DictReader(f))
    id_to_name = {
        str(row["accepted_host_id"]).strip(): norm_name(row["accepted_host_name"])
        for row in name_rows
    }
    all_ids = {hid for ids in panel_ids.values() for hid in ids}
    if set(id_to_name) != all_ids:
        raise RuntimeError("WCVP accepted-name map differs from frozen host ID set")

    if correction.get("schema", "").endswith("v0.2") and len(all_ids) != int(ident["accepted_host_ids"]):
        raise RuntimeError("v0.2 accepted-host ID breadth drift")

    host_species: dict[str, set[str]] = {}
    insect_hosts: dict[str, list[str]] = {}
    for insect in sorted(panel_ids):
        names = [id_to_name[hid] for hid in panel_ids[insect]]
        insect_hosts[insect] = names
        for host in names:
            genus = host.split()[0]
            host_species.setdefault(genus, set()).add(host)

    age_young, age_old = [int(x) for x in rule["time"]["primary_window_cal_yr_bp"]]
    dataset_types = [api_dataset_type(x) for x in rule["dataset"]["dataset_types"]]
    jobs = [(g, d) for g in sorted(host_species) for d in dataset_types]

    fetched: list[dict[str, Any]] = []
    failures: list[dict[str, str]] = []
    checkpoint_dir=(args.checkpoint_dir if args.checkpoint_dir is not None
                    else args.output_dir / "checkpoints")
    checkpoint_dir.mkdir(parents=True,exist_ok=True)
    missing_jobs=[]
    for genus,dtype in jobs:
        completed=read_completed_checkpoint(
            checkpoint_dir,genus,dtype,age_young,age_old
        )
        if completed is None:
            missing_jobs.append((genus,dtype))
        else:
            fetched.append(completed)
    print(json.dumps({
        "neotoma_frozen_queries":len(jobs),
        "reused_completed_queries":len(fetched),
        "remaining_queries":len(missing_jobs)
    }),flush=True)
    with ThreadPoolExecutor(max_workers=max(1, int(args.workers))) as pool:
        future_map = {
            pool.submit(fetch_genus, genus, dtype, age_young, age_old): (genus, dtype)
            for genus, dtype in missing_jobs
        }
        for fut in as_completed(future_map):
            genus, dtype = future_map[fut]
            try:
                result=fut.result()
                write_completed_checkpoint(
                    checkpoint_dir,genus,dtype,age_young,age_old,result
                )
                fetched.append(result)
            except Exception as exc:  # technical incompleteness, never absence
                failures.append({"genus": genus, "dataset_type": dtype, "error": str(exc)})
            print(json.dumps({
                "completed_unique_queries":len(fetched),
                "failed_queries":len(failures),
                "frozen_queries":len(jobs)
            }),flush=True)

    by_genus: dict[str, set[int]] = {g: set() for g in host_species}
    by_species: dict[str, set[int]] = {
        h: set() for names in host_species.values() for h in names
    }
    snapshot_rows: list[dict[str, Any]] = []
    for item in fetched:
        genus = item["genus"]
        dtype = item["dataset_type"]
        for taxonname, siteid in item["pairs"]:
            taxonname = norm_name(taxonname)
            siteid = int(siteid)
            by_genus[genus].add(siteid)
            if taxonname in by_species:
                by_species[taxonname].add(siteid)
            snapshot_rows.append(
                {
                    "query_genus": genus,
                    "dataset_type": dtype,
                    "taxonname": taxonname,
                    "siteid": siteid,
                }
            )

    exact_threshold = int(rule["coverage_census"]["exact_species_support_site_threshold"])
    genus_threshold = int(rule["coverage_census"]["genus_support_site_threshold"])

    host_rows: list[dict[str, Any]] = []
    for genus in sorted(host_species):
        genus_sites = len(by_genus[genus])
        for host in sorted(host_species[genus]):
            exact_sites = len(by_species[host])
            host_rows.append(
                {
                    "accepted_host_name": host,
                    "genus": genus,
                    "exact_species_lgm_sites": exact_sites,
                    "genus_lgm_sites": genus_sites,
                    "exact_species_supported": exact_sites >= exact_threshold,
                    "genus_supported": genus_sites >= genus_threshold,
                }
            )

    host_lookup = {r["accepted_host_name"]: r for r in host_rows}
    insect_rows: list[dict[str, Any]] = []
    for insect in sorted(insect_hosts):
        hosts = insect_hosts[insect]
        exact_n = sum(bool(host_lookup[h]["exact_species_supported"]) for h in hosts)
        genus_n = sum(bool(host_lookup[h]["genus_supported"]) for h in hosts)
        insect_rows.append(
            {
                "species": insect,
                "n_hosts": len(hosts),
                "exact_supported_hosts": exact_n,
                "genus_supported_hosts": genus_n,
                "any_exact_species_support": exact_n > 0,
                "any_genus_support": genus_n > 0,
                "all_hosts_genus_supported": genus_n == len(hosts),
            }
        )

    args.output_dir.mkdir(parents=True, exist_ok=True)
    snapshot_path = args.output_dir / "neotoma_lgm_taxon_site_snapshot.csv"
    with snapshot_path.open("w", newline="", encoding="utf-8") as f:
        fields = ["query_genus", "dataset_type", "taxonname", "siteid"]
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        w.writerows(sorted(snapshot_rows, key=lambda r: tuple(str(r[k]) for k in fields)))

    hosts_path = args.output_dir / "neotoma_host_coverage.csv"
    with hosts_path.open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(host_rows[0]))
        w.writeheader()
        w.writerows(host_rows)

    insects_path = args.output_dir / "neotoma_insect_coverage.csv"
    with insects_path.open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(insect_rows[0]))
        w.writeheader()
        w.writerows(insect_rows)

    summary = {
        "schema": "ttf_genetic_historical_host_neotoma_validation_result_v0.1",
        "status": (
            "PASS_RESPONSE_BLIND_NEOTOMA_COVERAGE_CENSUS"
            if not failures
            else "TECHNICALLY_INCOMPLETE_NEOTOMA_COVERAGE_CENSUS"
        ),
        "rule_sha256": sha256_path(args.rule),
        "panel_git_bytes_sha256_descriptive": sha256_path(args.panel),
        "identity_correction_sha256": sha256_path(args.identity_correction),
        "wcvp_names_sha256": sha256_path(args.wcvp_names),
        "logical_identity_digest_sha256": logical_digest,
        "primary_window_cal_yr_bp": [age_young, age_old],
        "dataset_types_api": dataset_types,
        "query": {
            "unique_host_species": sum(len(v) for v in host_species.values()),
            "unique_host_genera": len(host_species),
            "genus_x_dataset_queries": len(jobs),
            "completed_queries": len(fetched),
            "failed_queries": len(failures),
            "failures": failures,
        },
        "host_coverage": {
            "hosts": len(host_rows),
            "exact_species_supported": sum(bool(r["exact_species_supported"]) for r in host_rows),
            "genus_supported": sum(bool(r["genus_supported"]) for r in host_rows),
            "any_direct_lgm_site": sum(int(r["exact_species_lgm_sites"]) > 0 for r in host_rows),
            "any_genus_lgm_site": sum(int(r["genus_lgm_sites"]) > 0 for r in host_rows),
        },
        "insect_coverage": {
            "insects": len(insect_rows),
            "any_exact_species_support": sum(bool(r["any_exact_species_support"]) for r in insect_rows),
            "any_genus_support": sum(bool(r["any_genus_support"]) for r in insect_rows),
            "all_hosts_genus_supported": sum(bool(r["all_hosts_genus_supported"]) for r in insect_rows),
        },
        "outputs_sha256": {
            "snapshot": sha256_path(snapshot_path),
            "hosts": sha256_path(hosts_path),
            "insects": sha256_path(insects_path),
        },
        "decision_authority": False,
        "primary_panel_changed": False,
        "response_firewall": {
            "fresh_sequence_identity_opened": False,
            "fresh_pairwise_genetic_distances_opened": False,
            "fresh_post_ibd_turnover_computed": False,
            "fresh_beta_host_LGM_computed": False,
            "primary_palaeo_host_resistance_computed": False,
        },
    }
    summary_path = args.output_dir / "neotoma_validation_summary.json"
    summary_path.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    print(json.dumps(summary, sort_keys=True))
    return 0 if not failures else 2


if __name__ == "__main__":
    raise SystemExit(main())
