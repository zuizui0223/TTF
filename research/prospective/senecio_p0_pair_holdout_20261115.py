"""Exploratory P0 locality-pair phenotype score. Not a TTF species transfer claim.

Plan fixed in SENECIO_P0_PREDECLARED_LOCALITY_HOLDOUT_20261115.md
before reading phenotype values. Source SHA pinned; fail closed.
"""
from __future__ import annotations
import argparse
import csv
import hashlib
import io
import json
import math
import urllib.request
from pathlib import Path

COMMIT = "33222710aca75f79153e4118b409e7d77787d045"
REPO = "OB-lab/James_et_al._2021-Evolution"
PATH = "phenotype/data_files/traitMeansPairs.txt"
SHA = "400f73bf21fefb6fcb775396f256c3447c755595"
TRAITS = ("VH", "WW", "MSA", "MSD", "SBA", "AREA", "ELONG", "DISS", "CIRC")
ECOTYPES = ("Dune", "Headland")


def git_blob_sha(raw: bytes) -> str:
    return hashlib.sha1(f"blob {len(raw)}\0".encode() + raw).hexdigest()


def parse_pair_means(raw: bytes) -> dict:
    if git_blob_sha(raw) != SHA:
        raise ValueError("source SHA mismatch: fail closed")
    reader = csv.DictReader(io.StringIO(raw.decode("utf-8-sig")), delimiter="\t")
    fields = reader.fieldnames or []
    needed = {"Pair", "Ecotype"} | set(TRAITS)
    if len(fields) != len(set(fields)) or not needed <= set(fields):
        raise ValueError("missing or duplicate declared headers")
    parsed = {}
    total = 0
    for row in reader:
        total += 1
        pair, eco = (row["Pair"] or "").strip(), (row["Ecotype"] or "").strip()
        if not pair or eco not in ECOTYPES or (pair, eco) in parsed:
            raise ValueError("missing, invalid, or duplicated pair/ecotype ID")
        vals = []
        for name in TRAITS:
            number = float(row[name])
            if not math.isfinite(number):
                raise ValueError(f"non-finite {name}")
            vals.append(number)
        parsed[(pair, eco)] = tuple(vals)
    pairs = sorted({p for p, _ in parsed})
    if len(pairs) != 8 or total != 16 or any((p, e) not in parsed for p in pairs for e in ECOTYPES):
        raise ValueError("unexpected pair count or missing matching ecotype")
    return {p: tuple(parsed[(p, "Headland")][j] - parsed[(p, "Dune")][j] for j in range(len(TRAITS))) for p in pairs}


def score(raw: bytes) -> dict:
    changes = parse_pair_means(raw)
    pairs = sorted(changes)
    table = {}
    for j, trait in enumerate(TRAITS):
        vals = [changes[p][j] for p in pairs]
        preds = [sum(vals[k] for k in range(len(pairs)) if k != i) / (len(pairs) - 1) for i in range(len(pairs))]
        sse_zero = sum(y*y for y in vals)
        sse_test = sum((y-yhat)**2 for y, yhat in zip(vals, preds))
        agree = sum((y > 0 and yhat > 0) or (y < 0 and yhat < 0) for y, yhat in zip(vals, preds))
        ties = sum(y == 0 or yhat == 0 for y, yhat in zip(vals, preds))
        table[trait] = {
            "heldout_sse": sse_test,
            "zero_difference_sse": sse_zero,
            "heldout_skill": None if sse_zero == 0 else 1 - sse_test / sse_zero,
            "sign_agreements": agree,
            "sign_abstentions": ties,
            "sign_evaluable_pairs": len(pairs) - ties,
        }
    nonnull = sorted(v["heldout_skill"] for v in table.values() if v["heldout_skill"] is not None)
    median = (nonnull[(len(nonnull)-1)//2]+nonnull[len(nonnull)//2])/2 if nonnull else None
    return {
        "analysis": "exploratory_locality_pair_holdout_not_species_TTF",
        "source_repo": REPO, "source_commit": COMMIT, "source_path": PATH,
        "source_git_blob": SHA,
        "pair_ids": pairs,
        "n_independent_evaluation_units": len(pairs),
        "n_phenotype_traits": len(TRAITS),
        "metric": "leave_one_pair_out_mean_contrast_vs_zero_difference",
        "contrasts": "Headland-minus-Dune on published standardized trait means",
        "trait_results": table,
        "summary_descriptive": {
            "traits_with_positive_heldout_skill": sum(v["heldout_skill"] is not None and v["heldout_skill"] > 0 for v in table.values()),
            "median_trait_heldout_skill": median,
        },
        "claims_not_evaluated": [
            "cross_species_transfer", "genomic_mechanism", "local_adaptation",
            "reproductive_isolation", "historical_ancestry_additive_value",
        ],
    }


def source_bytes() -> bytes:
    url = f"https://raw.githubusercontent.com/{REPO}/{COMMIT}/{PATH}"
    with urllib.request.urlopen(url, timeout=30) as response:
        return response.read()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path)
    parser.add_argument("--output", type=Path)
    a = parser.parse_args()
    result = score(a.input.read_bytes() if a.input else source_bytes())
    body = json.dumps(result, ensure_ascii=False, sort_keys=True, indent=2) + "\n"
    if a.output:
        a.output.write_text(body)
    else:
        print(body)


if __name__ == "__main__":
    main()
