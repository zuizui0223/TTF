"""Source-only Senecio VCF header/sample-ID gate; no variant/genotype calls read.

No scientific response opening. Fetch pinned author tgz, verify Git blob SHA,
open header through #CHROM only, infer population keys conservatively.
"""
from __future__ import annotations
import argparse
import collections
import hashlib
import io
import json
import re
import tarfile
import urllib.request
from pathlib import Path

REPO = "OB-lab/James_et_al._2021-Evolution"
COMMIT = "33222710aca75f79153e4118b409e7d77787d045"
FILE = "genotype/vcf_file/ESC_rel_50pp_80md_HWE_pairs_MAF0.05.vcf.tgz"
GIT_BLOB_SHA = "fe1a452b183e0690e013cef57e59bedb13c7e072"
EXPECTED_POPULATIONS = (
    "D00", "D01", "D02", "D03", "D04", "D05", "D12", "D14", "D32",
    "H00", "H01", "H02", "H04", "H05", "H06", "H12", "H14", "H15"
)
MORPHOLOGY_PAIR_POPULATIONS = (
    ("D00", "H00"), ("D01", "H01"), ("D02", "H04"), ("D03", "H02"),
    ("D04", "H05"), ("D05", "H06"), ("D12", "H14"), ("D32", "H12"),
)
# Intentionally only a prefix. Never parse variant records.
SAMPLE_POPULATION = re.compile(r"^([DH][0-9]{2})(?:$|[^0-9])")


def git_blob_sha(raw: bytes) -> str:
    return hashlib.sha1(f"blob {len(raw)}\0".encode() + raw).hexdigest()


def read_vcf_header_only(fileobj) -> tuple[int, list[str]]:
    """Exit *immediately* after #CHROM, before first variant row."""
    meta_line_count = 0
    for index, raw in enumerate(fileobj, 1):
        if index > 50_000:
            raise ValueError("VCF #CHROM header not found within cap")
        if len(raw) > 5_000_000:
            raise ValueError("oversize VCF header line")
        if raw.startswith(b"##"):
            meta_line_count += 1
            continue
        if raw.startswith(b"#CHROM\t"):
            cols = raw.decode("utf-8", errors="strict").rstrip("\r\n").split("\t")
            if cols[:9] != ["#CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT"]:
                raise ValueError("unexpected VCF field layout")
            names = cols[9:]
            if not names or any(not x for x in names) or len(set(names)) != len(names):
                raise ValueError("missing or duplicate VCF sample IDs")
            return meta_line_count, names
        raise ValueError("unexpected non-header line before #CHROM")
    raise ValueError("VCF #CHROM not found")


def audit_archive(raw: bytes) -> dict:
    if git_blob_sha(raw) != GIT_BLOB_SHA:
        raise ValueError("source Git blob SHA mismatch")
    if len(raw) > 8_000_000:
        raise ValueError("unexpected compressed file size")
    with tarfile.open(fileobj=io.BytesIO(raw), mode="r:gz") as tf:
        members = [x for x in tf.getmembers() if x.isfile()]
        # macOS AppleDouble files are resource forks, not VCF content.
        vcfs = [x for x in members if x.name.endswith(".vcf")
                and not x.name.rsplit("/", 1)[-1].startswith("._")]
        if len(vcfs) != 1:
            names = sorted(x.name for x in members)[:16]
            raise ValueError(f"expected exactly one VCF among archived files; VCF count={len(vcfs)}, names={names}")
        member = vcfs[0]
        if member.size > 100_000_000 or member.size < 100:
            raise ValueError("unsafe or empty VCF decompressed length")
        with tf.extractfile(member) as stream:
            if stream is None:
                raise ValueError("missing VCF stream")
            n_meta, names = read_vcf_header_only(stream)
    counts = collections.Counter()
    malformed = 0
    for name in names:
        m = SAMPLE_POPULATION.match(name)
        if m and m.group(1) in EXPECTED_POPULATIONS:
            counts[m.group(1)] += 1
        else:
            malformed += 1
    paired = []
    for d,h in MORPHOLOGY_PAIR_POPULATIONS:
        paired.append({"morphology_pair": f"{d}-{h}",
                       "dune_sample_count": counts[d], "headland_sample_count": counts[h],
                       "both_population_prefixes_detected": bool(counts[d] and counts[h])})
    return {
        "status": "SOURCE_ONLY_VCF_HEADER_NO_GENOTYPE_OR_VARIANT_ROWS",
        "source_repo": REPO, "source_commit": COMMIT, "source_path": FILE,
        "source_git_blob": GIT_BLOB_SHA,
        "compressed_bytes": len(raw), "archive_member_name": member.name,
        "archive_member_declared_bytes": member.size, "vcf_meta_header_lines": n_meta,
        "vcf_sample_count": len(names),
        "population_prefix_pattern": SAMPLE_POPULATION.pattern,
        "population_sample_counts": dict(sorted(counts.items())),
        "unparsed_sample_count": malformed,
        "expected_population_prefixes_missing": sorted(set(EXPECTED_POPULATIONS) - set(counts)),
        "morphology_pair_population_coverage": paired,
        "morphology_pairs_with_both_population_prefixes": sum(p["both_population_prefixes_detected"] for p in paired),
        "phenotype_genotype_same_individuals": "NOT_VERIFIED",
        "snp_variant_values_read": False, "genotypes_read": False,
        "history_or_transfer_effect_estimated": False,
        "inference_gate": "NO_GENETIC_PREDICTION_AUTHORIZED",
    }


def download() -> bytes:
    u = f"https://raw.githubusercontent.com/{REPO}/{COMMIT}/{FILE}"
    with urllib.request.urlopen(u, timeout=45) as response:
        return response.read(8_000_001)


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--input", type=Path)
    p.add_argument("--output", type=Path)
    args = p.parse_args()
    result = audit_archive(args.input.read_bytes() if args.input else download())
    body = json.dumps(result, ensure_ascii=False, indent=2, sort_keys=True) + "\n"
    if args.output: args.output.write_text(body)
    else: print(body)


if __name__ == "__main__":
    main()
