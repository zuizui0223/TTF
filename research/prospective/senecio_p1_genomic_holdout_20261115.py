"""Exploratory P1: train-only SNP sign selection, heldout Senecio locality AUC.

Plan was committed before this code or first real VCF genotype opening:
SENECIO_P1_PRESPECIFIED_GENOMIC_HOLDOUT_20261115.md.
No full-data candidate/outlier lists or source phenotypes are used.
"""
from __future__ import annotations
import argparse
import collections
import hashlib
import io
import json
import random
import re
import statistics
import tarfile
import urllib.request
from pathlib import Path

SOURCE_COMMIT = "33222710aca75f79153e4118b409e7d77787d045"
SOURCE_REPO = "OB-lab/James_et_al._2021-Evolution"
SOURCE_PATH = "genotype/vcf_file/ESC_rel_50pp_80md_HWE_pairs_MAF0.05.vcf.tgz"
SOURCE_BLOB = "fe1a452b183e0690e013cef57e59bedb13c7e072"
PAIRS = [
    ("D00", "H00"), ("D01", "H01"), ("D02", "H04"), ("D03", "H02"),
    ("D04", "H05"), ("D05", "H06"), ("D12", "H14"), ("D32", "H12"),
]
EXTRA_PAIR = ("D14", "H15")
K = 100
MIN_POP_CALLS = 20
MIN_INDIVIDUAL_MARKERS = 70
REQUIRED_RECORDS = 9269
REQUIRED_SAMPLES = 1009
MISSING = 255
PATTERN = re.compile(r"^([DH][0-9]{2})(?:$|[^0-9])")


def git_blob_sha(data: bytes) -> str:
    return hashlib.sha1(f"blob {len(data)}\0".encode() + data).hexdigest()


def parse_gt(text: str) -> int:
    gt = text.split(":", 1)[0]
    sep = "|" if "|" in gt else "/" if "/" in gt else None
    if sep is None:
        return MISSING
    alleles = gt.split(sep)
    if len(alleles) != 2 or "." in alleles:
        return MISSING
    if not all(a in ("0", "1") for a in alleles):
        raise ValueError("GT allele unexpectedly exceeds one (not biallelic diploid)")
    return int(alleles[0]) + int(alleles[1])


def population_index(sample_names: list[str]) -> tuple[list[int], dict[str, list[int]], dict[str, int]]:
    if len(sample_names) != REQUIRED_SAMPLES or len(sample_names) != len(set(sample_names)):
        raise ValueError("unexpected original VCF sample header shape")
    expected = {p for pair in PAIRS + [EXTRA_PAIR] for p in pair}
    all_counts = collections.Counter()
    selected = []
    pop_inds = collections.defaultdict(list)
    for column, name in enumerate(sample_names):
        m = PATTERN.match(name)
        if m is None or m.group(1) not in expected:
            raise ValueError("unrecognized individual population prefix")
        pop = m.group(1)
        all_counts[pop] += 1
        if pop not in EXTRA_PAIR:
            pop_inds[pop].append(len(selected))
            selected.append(column)
    if len(all_counts) != 18 or len(pop_inds) != 16:
        raise ValueError("missing original population prefix")
    if len(selected) != 986:
        raise ValueError("morphology-linked genotype sample count differs from source-only audit")
    return selected, dict(pop_inds), dict(sorted(all_counts.items()))


def read_original_vcf(raw: bytes) -> tuple[list[dict], dict[str,list[int]], dict]:
    if git_blob_sha(raw) != SOURCE_BLOB:
        raise ValueError("source VCF tar Git blob hash mismatch")
    with tarfile.open(fileobj=io.BytesIO(raw), mode="r:gz") as tf:
        vcfs = [m for m in tf.getmembers() if m.isfile()
                and m.name.endswith(".vcf")
                and not m.name.rsplit("/",1)[-1].startswith("._")]
        if len(vcfs) != 1 or not 100 < vcfs[0].size < 100_000_000:
            raise ValueError("unique real VCF member is missing or unsafe")
        stream = tf.extractfile(vcfs[0])
        if stream is None:
            raise ValueError("VCF stream missing")
        sample_names = None
        selected = pop_idx = pop_counts = None
        records = []
        record_count = 0
        kept = 0
        unique_positions = set()
        for rawline in stream:
            if rawline.startswith(b"##"): continue
            if rawline.startswith(b"#CHROM\t"):
                fields = rawline.decode("utf-8").rstrip("\n\r").split("\t")
                if fields[:9] != ["#CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT"]:
                    raise ValueError("VCF record header schema mismatch")
                sample_names = fields[9:]
                selected, pop_idx, pop_counts = population_index(sample_names)
                continue
            if sample_names is None:
                raise ValueError("VCF variant line before #CHROM")
            record_count += 1
            if record_count > 12000: raise ValueError("more variants than expected VCF source")
            fields = rawline.decode("utf-8").rstrip("\n\r").split("\t")
            if len(fields) != 9+len(sample_names):
                raise ValueError("variant sample-column length mismatch")
            chrom, pos, _id, ref, alt = fields[:5]
            key = (chrom,pos)
            if key in unique_positions:
                raise ValueError("duplicate CHROM POS variant")
            unique_positions.add(key)
            if len(ref) != 1 or len(alt) != 1 or ref not in "ACGT" or alt not in "ACGT" or ref == alt:
                continue
            calls = bytearray(parse_gt(fields[9+i]) for i in selected)
            per_pop_freq, per_pop_n = {}, {}
            for pop, indices in pop_idx.items():
                good = [calls[i] for i in indices if calls[i] != MISSING]
                per_pop_n[pop] = len(good)
                per_pop_freq[pop] = (sum(good)/(2*len(good))) if good else None
            records.append({"chrom":chrom, "pos":int(pos),
                            "calls":calls, "freq":per_pop_freq, "n":per_pop_n})
            kept += 1
        if record_count != REQUIRED_RECORDS:
            raise ValueError(f"expected {REQUIRED_RECORDS} original VCF variant records, found {record_count}")
        if not records:
            raise ValueError("no valid biallelic SNP")
    return records, pop_idx, {
        "source_blob": SOURCE_BLOB, "vcf_header_sample_count":len(sample_names),
        "population_counts":pop_counts, "matched_sample_count":len(selected),
        "source_variant_records":record_count, "biallelic_snv_records":kept,
        "distinct_snv_contigs":len({r["chrom"] for r in records}),
    }


def train_markers(records: list[dict], omitted: int) -> list[tuple]:
    training_pairs = [p for i,p in enumerate(PAIRS) if i != omitted]
    best_by_contig = {}
    for locus in records:
        if not all(locus["n"][pop] >= MIN_POP_CALLS for pair in training_pairs for pop in pair):
            continue
        delta = [(locus["freq"][h]-locus["freq"][d]) for d,h in training_pairs]
        mean_delta = sum(delta)/len(delta)
        if mean_delta == 0.0: continue
        f0 = sum(locus["freq"][pop] for pair in training_pairs for pop in pair) / 14
        rank = (-abs(mean_delta), locus["chrom"], locus["pos"])
        candidate = (rank, locus, 1 if mean_delta > 0 else -1, f0)
        key = locus["chrom"]
        if key not in best_by_contig or rank < best_by_contig[key][0]:
            best_by_contig[key] = candidate
    selected = sorted(best_by_contig.values(), key=lambda x:x[0])[:K]
    return selected


def auc_of(scores_a: list[float], scores_b: list[float]) -> float:
    """Headland scores A vs Dune scores B, pairwise Mann-Whitney with ties."""
    if not scores_a or not scores_b:
        raise ValueError("AUC has empty ecotype group")
    win = 0.
    for x in scores_a:
        for y in scores_b:
            win += 1.0 if x > y else 0.5 if x == y else 0.
    return win / (len(scores_a)*len(scores_b))


def score_one_fold(selected: list[tuple], pop_idx: dict[str,list[int]], hold: tuple[str,str]) -> dict:
    if len(selected) != K:
        return {"status":"NOT_EVALUABLE_INSUFFICIENT_CONTIGS",
                "contig_representatives":len(selected)}
    d,h=hold
    results={}
    for pop in (d,h):
        scores=[]
        n_calls=[]
        for idx in pop_idx[pop]:
            z=0.
            count=0
            for _, locus, direction, f0 in selected:
                value=locus["calls"][idx]
                if value == MISSING:
                    # Training-only reference mean imputation gives zero contribution.
                    continue
                count += 1
                z += direction*(value/2-f0)
            n_calls.append(count)
            if count >= MIN_INDIVIDUAL_MARKERS:
                scores.append(z/K)
        results[pop]={"scores":scores, "n_total":len(pop_idx[pop]),
                      "n_evaluable":len(scores),
                      "mean_marker_call_rate":sum(n_calls)/(K*len(n_calls))}
    if min(results[d]["n_evaluable"],results[h]["n_evaluable"]) < MIN_POP_CALLS:
        return {"status":"NOT_EVALUABLE_TEST_CALLS", "pair":f"{d}-{h}",
                "test_population_called":{pop:{k:v for k,v in results[pop].items() if k!="scores"} for pop in (d,h)}}
    aligned=ties=eligible=0
    for _,locus,direction,f0 in selected:
        if min(locus["n"][d],locus["n"][h]) < MIN_POP_CALLS:
            continue
        val=locus["freq"][h]-locus["freq"][d]
        if val==0:
            ties+=1
            continue
        eligible+=1
        aligned+=int((val>0)==(direction>0))
    return {
        "status":"EVALUABLE_LOCALITY_AUC",
        "pair":f"{d}-{h}",
        "auc":auc_of(results[h]["scores"], results[d]["scores"]),
        "n_headland_test":results[h]["n_evaluable"],
        "n_dune_test":results[d]["n_evaluable"],
        "n_headland_original":results[h]["n_total"],
        "n_dune_original":results[d]["n_total"],
        "headland_marker_call_rate":results[h]["mean_marker_call_rate"],
        "dune_marker_call_rate":results[d]["mean_marker_call_rate"],
        "signed_snp_direction_agreement_count":aligned,
        "signed_snp_direction_evaluable_count":eligible,
        "signed_snp_direction_ties":ties,
        "signed_snp_direction_fraction":aligned/eligible if eligible else None,
        "distinct_training_selected_contigs":len({locus["chrom"] for _,locus,_,_ in selected}),
    }


def bootstrap_pair_auc(auc:list[float], n:int=2000, seed:int=20261115) -> list[float]:
    rand=random.Random(seed)
    boot=sorted(sum(auc[rand.randrange(len(auc))] for _ in auc)/len(auc) for i in range(n))
    return [boot[int(.025*(n-1))],boot[int(.975*(n-1))]]


def analyze(raw: bytes) -> dict:
    records, pop_idx, source = read_original_vcf(raw)
    folds=[]
    for hold_index, pair in enumerate(PAIRS):
        chosen=train_markers(records,hold_index)
        folds.append(score_one_fold(chosen,pop_idx,pair))
    failed=[x for x in folds if x["status"] != "EVALUABLE_LOCALITY_AUC"]
    if failed:
        return {"status":"NOT_EVALUABLE_PREDECLARED_GENOMIC_GATE",
                "source":source, "folds":folds, "failure_pairs":[x.get("pair") for x in failed],
                "no_posthoc_rescue":True}
    auc=[x["auc"] for x in folds]
    return {
        "status":"EXPLORATORY_P1_GENOMIC_LOCALITY_TRANSFER_REPLAY",
        "source_repo":SOURCE_REPO, "source_commit":SOURCE_COMMIT,
        "source_path":SOURCE_PATH, "source_git_blob":SOURCE_BLOB,
        "source":source,
        "heldout_unit":"one Dune-vs-Headland population pair",
        "heldout_pair_count":len(PAIRS),
        "train_only_marker_count_per_fold":K,
        "marker_pruning":"one best training SNP per VCF contig",
        "folds":folds,
        "summary":{
            "pair_equal_mean_auc":sum(auc)/len(auc),
            "min_auc":min(auc),"max_auc":max(auc),
            "pair_equal_median_auc":statistics.median(auc),
            "auc_pair_cluster_bootstrap_2000_percentile":bootstrap_pair_auc(auc),
            "pairs_auc_above_0_5":sum(x>0.5 for x in auc),
        },
        "source_asertainment_limitation": "original source VCF SNP filters used all 18 populations",
        "not_claimed":["species_transfer","causal_genetic_loci","adaptive_selection",
                       "reproductive_isolation","CREST_history_effect","independent_species_replication"],
        "no_raw_genotypes_exported":True,
        "analysis_branch_is_separate_from_frozen_TTF_science":True,
    }


def download_original() -> bytes:
    u=f"https://raw.githubusercontent.com/{SOURCE_REPO}/{SOURCE_COMMIT}/{SOURCE_PATH}"
    with urllib.request.urlopen(u, timeout=40) as resp:
        return resp.read(8_000_001)


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--input",type=Path)
    p.add_argument("--output",type=Path)
    args=p.parse_args()
    result=analyze(args.input.read_bytes() if args.input else download_original())
    body=json.dumps(result, indent=2, sort_keys=True, ensure_ascii=False) + "\n"
    if args.output:
        args.output.write_text(body)
    else:
        print(body)


if __name__=="__main__":
    main()
