# Penstemon original Table S1 visual amendment, November 18, 2026

**Status: one structural source discrepancy RESOLVED, biological phenotype availability HOLD remains.**

The original publisher-authoritative 21-page Europe PMC supplement `media-1.pdf` (SHA256 `3682720dad4f0621bebef3be7014aa667067c5f0dce9160ca8f4c4a758e7f943`) was rendered as original, unmodified page-11 (Table S1) and top-of-page-12 (Table S2) PNGs after the same-authority download retry succeeded. The independent original source run is [Actions #44300695234](https://github.com/zuizui0223/TTF/actions/runs/44300695234), artifact ID **14426680096**, artifact ZIP SHA256 **0799536c8374336c3f9cdbdbe1659ac4b7b5f78d547e382300e4eaed14972728**.

## One documented original formatting exception

The Table S1 visual contains an unhyphenated voucher **`TTC026503`**, captioned **`P. cardinalis`**, marked **`Phenotyped=no`**, with the original collector `MJ` and herbarium `TTC`. All previous 51 correctly extracted rows used an uppercase alphabetic prefix separated from digits by a hyphen. Thus this single, exact visually authenticated source exception explains the missing **52nd voucher** and the former shortage of one *P. cardinalis* relative to NCBI's 2 BioSample taxon labels.

**Parser repair:** accept exactly `TTC026503` in addition to the original hyphenated voucher grammar. Do NOT permit arbitrary unhyphenated IDs or reclassify that row to phenotyped. The original strict 52-voucher/**34 phenotyped** admission gate is unchanged. Its new expected source result is **52 total vouchers, only 30 author `Phenotyped=yes` flags**, so this is NOT a full phenotype-source authorization.

## What is still unresolved

1. Original prose reports **34 phenotyped accessions**, while the visually verified Table S1 contains **30 `yes` rows** and 22 `no` rows after adding `TTC026503`. The article also reports four phenotyped accessions measured in the field for floral dimensions; **the relation to S1 `yes/no` must be proven rather than assumed**.
2. Original Table S2 starts at the **top of PDF page 12** (visual evidence). It is explicitly titled *Information for populations phenotyped for floral dimension traits in the field*. Thus the next safe check is a **source-ID-only** listing of S2 fields, matching them conservatively to S1 voucher IDs without accessing quantitative flower traits.
3. The original Table S1 visibly names voucher `CAW-190` as `P. strictiformis` while NCBI includes two `Penstemon comarrhenus` entries but only one literal match in the S1 text. Synonym/identification equivalence has **not** been independently verified and is not silently repaired.
4. The S1 voucher IDs are **not NCBI SAMN accessions**. The 52 GenBank BioSample identity mapping remains taxon-level until an authentic same-voucher mapping is found.
5. Table S3 in the same original supplement is phylogenetic ANOVA model results, **not a raw accession-by-trait table**. The successful source-only S3 check [Actions #43875710722](https://github.com/zuizui0223/TTF/actions/runs/43875710722) found zero voucher keys, so it cannot close the 34-genotype phenotype join.

## TTF-Q consequence

The four original hummingbird-syndrome origin groups and their bee sister taxa are well grounded in **52 sequenced BioSample species labels**. No source-matched four-fold TTF-P, historical-genetic increment, or reproductive isolation prediction has yet been computed. Source-only S1 and S2 measurement availability can be profiled, but **retain `PHENOTYPE_34_SAMN_JOIN_NOT_VERIFIED`**, not a biological null and not an authorized confirmatory out-of-origin prediction.

Scientific freeze on main remains intact.
