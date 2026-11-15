# Senecio P1: prespecified exploratory genomic-ecotype transfer across localities (2026-11-15)

**Registration point:** The P1 plan is committed before reading any VCF variant/genotype records. P0 phenotype outcomes were already opened on separate PR #157. P1 does not use phenotype outcomes for marker selection, thresholds, or locality choice; no frozen TTF-Q, COI, He, Fst or break-location test is reopened. This is NOT a confirmatory TTF application and not a cross-species transfer claim.

## Original source and comparison units

- Original James et al. (2021), Evolution, author repository https://github.com/OB-lab/James_et_al._2021-Evolution , commit 33222710aca75f79153e4118b409e7d77787d045.
- Original archive: genotype/vcf_file/ESC_rel_50pp_80md_HWE_pairs_MAF0.05.vcf.tgz ; Git blob SHA fe1a452b183e0690e013cef57e59bedb13c7e072.
- Verified **source-only** in PR #156, Actions #43691834787: 1,009 distinct VCF sample IDs, 18 populations, 986 individuals across 16 populations that form eight morphology-linked locality pairs. The 23 individuals in D14-H15 have no paired phenotype source and are excluded here.
- Eight independent heldout *locality-pair* units: D00-H00, D01-H01, D02-H04, D03-H02, D04-H05, D05-H06, D12-H14, D32-H12.
- Prespecified source validation: exact original Git blob SHA; 1,009 sample headers; 9,269 original VCF records; exact 18 population prefixes; no missing or duplicate sample IDs. Only biallelic A/C/G/T SNVs considered. Treat missing/partial diploid GT as missing, fail on unexpected higher allele number or dosage representation.

## Primary question and exact model

Can a *signed genomic contrast* learned in seven locality pairs discriminate Headland (positive ecotype) from Dune genotypes in the entirely withheld eighth locality pair?

For locus l in pair p, d[p,l] is alternate allele frequency in Headland minus that in Dune, estimated from observed diploid genotypes.

**Train-only feature and direction selection (eight leave-one-pair-out folds):**

1. For each fold, within each of the 14 training populations require at least 20 observed diploid genotype calls for a locus. No heldout allele frequency, Fst or genotype is allowed into feature ranking.
2. Compute m[l] = mean of the seven pair-equal signed allele-frequency differences d[p,l].
3. Among loci sharing the same VCF contig (CHROM), keep the largest absolute m[l]; lexicographic CHROM/POS tie-break. Rank these contig representatives by absolute m[l] and choose exactly K=100. If fewer than 100 contigs are eligible, stop with NOT_EVALUABLE_INSUFFICIENT_CONTIGS. K is immutable.
4. For each selected locus retain sign(m[l]) and f0[l], the equal-weighted mean of the 14 training population allele frequencies.
5. Each heldout individual receives a score z = mean over the 100 selected loci of sign(m[l]) times (alternate dosage/2 − f0[l]). For a missing genotype, use the training f0[l] so that its contribution is zero. No cutoff or scaling is learned from the heldout labels.
6. Test individuals must have at least 70 callable selected loci. If either heldout ecotype has fewer than 20 admitted individuals, stop without computing AUC for that fold.

**Primary statistic:** tie-corrected heldout AUC for Headland versus Dune based on z, reported individually for each of eight pairs. The null AUC=0.5 is only a descriptive reference, not a null-hypothesis p-value. Summarize the eight pair-equal mean AUC and descriptive percentile interval from 2,000 pair-cluster bootstrap resamples with seed 20261115. Small nonrandom n=8 limits generality.

**Secondary statistic:** among selected 100 contig-representative SNPs with at least 20 called individuals in each heldout ecotype and nonzero heldout signed frequency contrast, the fraction for which the heldout sign agrees with the seven-pair training sign. Report denominator and the selected-marker missingness by heldout ecotype. Contig pruning does not prove complete linkage independence. No pooled SNP-level significance test.

## Claim firewall

- No source-author all-pair outlier SNP list can be used for predictor choice: all marker ranking and sign estimation are recomputed within each seven-pair training fold.
- Source SNP ascertainment (MAF/missingness filtering across original 18 populations) itself used the wider source panel and is not train-only. Report as a limitation of interpreting out-of-locality transfer; a raw-sequence independent ascertainment test would be stronger.
- The eight comparisons are one species, not eight independently speciating taxa. Any positive AUC does not prove causal alleles, selection, gene function, historical ancestry effects, or reproductive isolation.
- Shared ancestry, regional structure, uneven genotype missingness, reference orientation, and sampling overlap can mimic cross-locality genomic predictability. Report missingness and do not equate AUC with mechanism homology.
- Pair-level P0 morphological transfer remains separate; genotype individual IDs have not been verified as the same as phenotype plants. The two are comparable at the population-pair level only.
- No after-the-fact feature count, missingness threshold, model family, SNP orientation, selection strategy or subgroup can rescue the registered P1 result. If a gate fails, report NOT_EVALUABLE.
- No raw VCF/individual genotypes will be committed, only auditable aggregated JSON and code.

## Planned assets

An exact-source-pinned stdlib Python replay, synthetic leakage tests, and a one-shot GitHub Actions report. No TTF main branch edits.
