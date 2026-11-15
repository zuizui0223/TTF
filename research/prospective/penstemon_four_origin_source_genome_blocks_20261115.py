"""TTF-Q four independent origin BLOCK holdout from verified source BioSamples.

Only authentic 52 BioSample/taxon labels from frozen 2026 ENA+NCBI receipt.
No flowers, SNP genotypes, inherited-state estimates or reproductive isolation
effects have been opened; output contains only source support and worst-case
combinatorial phenotype-access bounds. Protocol preregistered in
PENSTEMON_FOUR_ORIGIN_SOURCE_ONLY_BLOCKS_20261115.md.
"""
from __future__ import annotations
import collections
import hashlib
import json
from pathlib import Path

SOURCE=Path(__file__).with_name("penstemon_habroanthus_52_ncbi_biosample_taxon_receipt_20261115.json")
BLOB_SHA="acdf2997e4d338fcb879d0323dc4c6a54f4fe413"
REPORTED_PHE_NO=34
REPORTED_GENOMES=52
REPORT_TAXA=21
BLOCKS={
    "B1_barbatus_virgatus":{
        "derived":"Penstemon barbatus",
        "sisters":["Penstemon virgatus"],
    },
    "B2_cardinalis_glaber":{
        "derived":"Penstemon cardinalis",
        "sisters":["Penstemon glaber"],
    },
    "B3_eatonii_laevis":{
        "derived":"Penstemon eatonii",
        "sisters":["Penstemon laevis"],
    },
    "B4_labrosus_deaveri_putus":{
        "derived":"Penstemon labrosus",
        "sisters":["Penstemon deaveri","Penstemon putus"],
    },
}


def source_blob_sha(b:bytes)->str:
    return hashlib.sha1(f"blob {len(b)}\0".encode()+b).hexdigest()


def specimen_geometry(d:dict)->dict:
    if d.get("status")!="PENSTEMON_HABROANTHUS_52_ORIGINAL_BIOSAMPLE_SOURCE_TAXON_GATE":
        raise ValueError("unexpected source BioSample receipt status")
    s=d.get("source_matched_taxonomy") or {}
    if s.get("status")!="ALL_52_ORIGINAL_BIOSAMPLES_TAXON_SOURCE_MATCHED":
        raise ValueError("all original 52 BioSample taxon identities not authenticated")
    rows=s.get("source_biosample_taxon_crosswalk") or []
    ids=[r.get("biosample_accession") for r in rows]
    taxa=[r.get("scientific_name") for r in rows]
    if len(rows)!=REPORTED_GENOMES or len(set(ids))!=REPORTED_GENOMES:
        raise ValueError("source specimen ID count/uniqueness drift")
    if len(set(taxa))!=REPORT_TAXA or any(not x for x in ids+taxa):
        raise ValueError("source taxon labels missing or taxonomy count drift")
    if collections.Counter(taxa)!=s.get("source_taxon_sample_counts"):
        raise ValueError("source taxon count and source accession rows disagree")
    allmembers=set()
    n_total=0
    blocks={}
    for key,meta in BLOCKS.items():
        members=[meta["derived"]]+list(meta["sisters"])
        if len(members)!=len(set(members)) or allmembers.intersection(members):
            raise ValueError("overlapping/duplicated biological-origin taxon groups")
        allmembers.update(members)
        counts={name:taxa.count(name) for name in members}
        if any(n<1 for n in counts.values()):
            raise ValueError(f"source-defined sister comparator lacking original genome: {key}")
        block_count=sum(counts.values())
        n_total += block_count
        # Sharp bounds when 34 phenotyped accessions are an unknown
        # subset of 52 genome specimens, not a random subset.
        blocks[key]={
            "origin_derived_species":meta["derived"],
            "origin_source_derived_genome_samples":counts[meta["derived"]],
            "source_sister_taxa":{t:counts[t] for t in meta["sisters"]},
            "genome_sample_count":block_count,
            "true_phenotyped_genome_sample_count":"NOT_SOURCE_MATCHED",
            "possible_phenotyped_count_bounds_from_34_of_52":[
                max(0, block_count-(REPORTED_GENOMES-REPORTED_PHE_NO)),
                min(block_count,REPORTED_PHE_NO)
            ],
            "source_genomic_comparison_available":True,
        }
    if n_total>len(rows):
        raise ValueError("origin blocks expand sample total beyond original 52")
    folds=[]
    for heldout,detail in blocks.items():
        training={k:v for k,v in blocks.items() if k!=heldout}
        if len(training)!=3:
            raise ValueError("exact 3-train/1-evaluation origin folds not formed")
        folds.append({
            "heldout_origin_block":heldout,
            "training_origin_block_count":len(training),
            "training_genome_samples":sum(v["genome_sample_count"] for v in training.values()),
            "heldout_genome_samples":detail["genome_sample_count"],
            "training_derived_genome_samples":sum(v["origin_source_derived_genome_samples"] for v in training.values()),
            "training_bee_comparator_genome_samples":sum(sum(v["source_sister_taxa"].values()) for v in training.values()),
            "heldout_derived_genome_samples":detail["origin_source_derived_genome_samples"],
            "heldout_bee_comparator_genome_samples":sum(detail["source_sister_taxa"].values()),
            "heldout_phenotype_measurements":"NOT_VERIFIED",
        })
    outside_count=REPORTED_GENOMES-n_total
    return {
        "status":"FOUR_ORIGIN_GENOMIC_SOURCE_COVERAGE_ONLY_PHENOTYPE_JOIN_HOLD",
        "source_biosample_git_blob":BLOB_SHA,
        "bioProject":"PRJNA1479764",
        "source_total_genome_samples":REPORTED_GENOMES,
        "source_total_taxon_labels":REPORT_TAXA,
        "paper_reported_pheno_accessions":REPORTED_PHE_NO,
        "author_derived_origin_blocks":len(BLOCKS),
        "strict_origin_heldout_folds":folds,
        "block_support":blocks,
        "genome_samples_in_four_origin_blocks":n_total,
        "genome_samples_outside_origin_blocks":outside_count,
        "phenotyped_samples_in_four_origin_blocks_possible_bounds":[
            max(0, REPORTED_PHE_NO-outside_count),
            min(REPORTED_PHE_NO,n_total),
        ],
        "phenotype_sample_ID_join":"NOT_VERIFIED",
        "genome_sample_taxon_accession_join":"VERIFIED_SOURCE",
        "unseen_evolutionary_origin_evaluation_unit_count":4,
        "no_flower_trait_snp_or_fitness_outcome_read":True,
        "no_TTF_cross_origin_result_claim":True,
    }


def main():
    raw=SOURCE.read_bytes()
    if source_blob_sha(raw)!=BLOB_SHA:
        raise ValueError("original source-BioSample receipt git blob drift")
    print(json.dumps(specimen_geometry(json.loads(raw)),sort_keys=True,indent=2))


if __name__=="__main__":main()
