# Senecio lautus: cross-layer locality transfer (P0 morphology versus P1 genome)

**Data date:** 2026-11-15; author original James et al. (2021), Evolution. **Evidence status:** authentic source-backed but EXPLORATORY, within one species, not a confirmatory TTF-Q result, not an independent field replication, not a test of causal genomic mechanism or speciation.

## Preregistered original genotype result

The P1 design (7 training locality pairs → 1 untouched locality) was committed before original genotype rows were read in [P1 prospective design](SENECIO_P1_PRESPECIFIED_GENOMIC_HOLDOUT_20261115.md). Synthetic no-leakage tests passed at Actions [#43693185386](https://github.com/zuizui0223/TTF/actions/runs/43693185386). The source-opening event was one-shot, SHA-locked and succeeded at [Actions #43693277926](https://github.com/zuizui0223/TTF/actions/runs/43693277926). Exact aggregate receipt: [senecio_p1_first_opening_aggregate_receipt_20261115.json](senecio_p1_first_opening_aggregate_receipt_20261115.json). Only the 16 morphology-linked populations were used, with 986 distinct VCF sample IDs; original 9,269 biallelic SNP records across 3,431 distinct contigs, selected 100 contig-representative SNPs from the seven training locality pairs in each fold. The morphology-free ninth genomic D14–H15 pair was excluded.

- **Primary result:** Pair-equal eight-fold held-out AUC **0.463809**; 2,000 pair-bootstrap percentile interval **0.180303–0.771421** (descriptive, not a generalizable 95% population confidence interval); median **0.353486**; **3/8** locality pairs above 0.5.
- Per-population marker call rate among heldout samples was approximately 94.5–97.7%, and all 986 original test samples passed the fixed >=70/100 marker callability gate. Thus a gross difference in overall marker availability does not by itself explain the extreme split; finer genotyping/batch structure remains unassessed.
- The 100 selected contig-SNP allele-frequency directions did not replicate coherently in every withheld pair: informative selected-SNP sign concordance fractions across pairs range approximately **0.405–0.568**. AUC may be near zero even when pairwise SNP sign concordance is near 0.5, because the effect magnitudes and covariance of genotype scores matter.

## Post-outcome descriptive comparison to previously opened phenotype P0

The separate P0 study [PR #157](https://github.com/zuizui0223/TTF/pull/157) had found 9/9 original morphology traits improving heldout locality-pair mean-difference prediction versus a fixed zero-difference model. We now take the **same eight locality pair names**, recompute each P0 heldout contrast sign for each of nine author-standardized traits and compare to the frozen P1 AUC. This is explicitly **post-outcome descriptive**, not a predeclared joint hypothesis test.

| Locality pair | P0 morphology sign agreement (of 9 traits) | P1 genotype individual classification AUC |
|---|---:|---:|
| D00–H00 | 8 | 0.956477 |
| D01–H01 | 6 | 0.983046 |
| D02–H04 | 8 | 0.973985 |
| D03–H02 | 9 | 0.007525 |
| D04–H05 | 8 | 0.082466 |
| D05–H06 | 8 | 0.361863 |
| D12–H14 | 9 | 0.000000 |
| D32–H12 | 9 | 0.345110 |

Five pairs (**D03–H02, D04–H05, D05–H06, D12–H14, D32–H12**) have 8–9/9 P0 morphological contrast directions agreeing with the leave-pair-out training trend, but the P1 signed genomic ecotype score has AUC < 0.5. Two comparisons have AUC <0.1 and one is exactly zero: the trained sign points opposite to the heldout ecotype labels rather than simply losing all information. The P0 and P1 score scales/estimands differ (traits' mean contrasts versus individual-genotype AUC), so their numerical magnitudes cannot be subtracted or pooled. [Machine-readable postoutcome cross-layer table](senecio_p0_p1_crosslayer_descriptive_20261115.json).

**Conservative biological interpretation:** repeatable dune-versus-headland morphological direction is compatible with heterogeneous or reversed signed genomic score portability across locality populations. This does not imply distinct causal alleles, different genetic developmental mechanisms, adaptive genomic novelty, different speciation stages, or known MLTR historical state modes. Geographic ancestry, neutral structure, original all-population SNP ascertainment, linkage across contigs, and metadata/collection biases remain competing explanations.

## What the result does and does not license

**Supported within the declared source:** The one fixed, seven-pair-trained, 100-contig SNP directional score is not consistently oriented in an entirely heldout eighth locality, despite a separate morphological contrast signal in the same set of populations. The result is a nontrivial observational mismatch between two **different predictive interfaces**.

**Not yet supported:** (1) an ancestral-history predictor that improves prediction; (2) a species-level out-of-domain TTF result (all observations are one Senecio species); (3) a molecular cause for any morphological ecotype trait; (4) a reproductive-isolation outcome; (5) a general CREST/MLTR theorem validated in nature.

**Next biological discriminator:** independently constructed neutral population ancestry/IBD and historical demographic structure, *without using the selected SNPs as both the explanatory history predictor and the response*, followed by an independent population/lineage validation. With only eight locality pairs, any additional ancestry-predictor fit remains exploratory and must be predeclared; external Timema or Costus species panels are still needed for actual species-holdout TTF.

## Source and integrity

- Original VCF SHA: fe1a452b183e0690e013cef57e59bedb13c7e072, author commit 33222710aca75f79153e4118b409e7d77787d045.
- Original phenotype pair-mean table Git blob: 400f73bf21fefb6fcb775396f256c3447c755595, corresponding to 16 population means.
- P1 original output was not subsequently retuned or re-opened. Authors' full-data candidate SNP outlier lists are not used for the signed predictor. However original VCF SNP discovery/MAF filters were performed before this cross-validation on the complete dataset, imposing ascertainment conditioning.
- Aggregate-only artifact from the one-shot GitHub workflow: ID 14153331627, SHA256 99acc1717fba2a19e406d3c2a09e7beb1b9ce1871e726f68448bad2c530363ec. Raw VCF and individual sample genotypes are not redistributed.
