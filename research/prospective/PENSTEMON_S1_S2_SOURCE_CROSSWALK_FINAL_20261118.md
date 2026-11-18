# Penstemon original source S1–S2: verified voucher availability and strict TTF-Q hold (2026-11-18)

**Result:** exact original source PDF voucher **52/52** recovered, not 51/52, via visually verified unhyphenated `TTC026503` (*P. cardinalis*, `Phenotyped=no`). However manuscript's **34 phenotyped** is still not identical to the **30 yes + 22 no** literal original Table S1 flags; this is unresolved study-source provenance, **not** a null effect.

## Official source and direct check

Stone et al. (2026) preprint DOI [10.64898/2026.07.21.739878](https://doi.org/10.64898/2026.07.21.739878), Europe PMC source supplement `media-1.pdf`, SHA256 `3682720dad4f0621bebef3be7014aa667067c5f0dce9160ca8f4c4a758e7f943`, authenticated full-page image [GitHub Actions #44300695234](https://github.com/zuizui0223/TTF/actions/runs/44300695234), artifact 14426680096. Original source-only S1/S2 ID reparse independently passed [Actions #44302939000](https://github.com/zuizui0223/TTF/actions/runs/44302939000), aggregate-only artifact 14427815767, ZIP SHA256 `1529c313593ff618d78131890c2fb973a5f281c994b2e815bf8d45de701f8d70`.

| Source support | Exact verified count |
|---|---:|
| NCBI original whole-genome BioSamples | 52 |
| Table S1 distinct source voucher IDs | **52** |
| Table S1 `Phenotyped=yes` | **30** |
| Table S1 `Phenotyped=no` | **22** |
| Manuscript phenotyped source total | **34 reported, not reconciled** |
| Original S2 independent field floral-dimension source IDs | **25** |
| S2 field IDs that **also occur** as S1 genomic specimen voucher IDs | **5** |
| S2 field IDs not in S1 | **20** |
| Original S1 taxon literal not in 21 source NCBI labels | **1, `CAW-190`** |
| Original S3 accession-by-trait keys | **0** (S3 is phylogenetic ANOVA results, not raw measurements) |

The **five exact S2→S1 crosswalks** are `CAW015`-style IDs normalized only by punctuation:
`CAW080↔CAW-080`, `CAW167↔CAW-167`, `CAW176↔CAW-176`, `CAW180↔CAW-180`, `CAW181↔CAW-181`. This is a **source voucher-key overlap**, not guaranteed that the same physical plant or the same sequenced individual was phenotyped. Twenty original S2 rows are **not** part of the 52-genome S1 set. S2 consists of field flower-dimensional measurement populations; it is not a hidden 4-sample correction to the 34-genome-individual denominator.

The visually verified `CAW-190` row has literal *P. strictiformis*, while the 52 NCBI source taxonomy includes two *P. comarrhenus*. **No synonym or accession identity is inferred**. Without an author-source synonym map and explicit voucher-to-SAMN mapping, a fully authenticated 52 genotype × 34 phenotype individual join is still **HOLD**.

## Four independently inferred pollination origins: S1 source `yes` coverage

| Original origin block | Derived-lineage S1 `yes` | Bee comparator(s) S1 `yes` | Both source phenotyped groups? |
|---|---:|---:|---|
| *P. barbatus* ↔ *P. virgatus* | 4 | 4 | Yes |
| *P. cardinalis* ↔ *P. glaber* | 1 | 1 | Yes |
| *P. eatonii* ↔ *P. laevis* | 2 | 1 | Yes |
| *P. labrosus* ↔ *P. deaveri* + *P. putus* | 2 | 3 + 2 | Yes |

All four **paper-defined** origins are represented by at least one derived and one sister lineage with the original S1 `Phenotyped=yes` source flag; this is source-level availability, **not** a valid species-heldout genomic/trait transport score. Realized sample count per independent origin is very small, and there are only **four** evolutionary units regardless of the 52 sampled genomes. Same-voucher/trait panel and validated original phenotyped 34 remain essential.

## TTF/CREST admission decision

- `GENOME_SOURCE_52_PASSED`, `S1_VOUCHER_52_PASSED`, `S1_YES_30_PASSED_AS_LITERAL_NOT_AS_34`, `S2_SOURCE_25_5OVERLAP_PASSED`.
- `TAXON_LITERAL_IDENTITY_CA_W190_UNRESOLVED`, `PHENOTYPED_34_CONCORDANCE_UNRESOLVED`, `S1_VOUCHER_TO_SAMN_UNVERIFIED`.
- **TTF-P across independent evolutionary origins = NOT_EVALUABLE** until strict source phenotype and lineage IDs are fixed. Original study qualitative phenotypic convergences should not be presented as our new biological discovery.

[Reproducible source code](penstemon_original_s1_52_s2_field_key_crosswalk_20261118.py) · [frozen original source ID receipt](penstemon_2026_original_S1_52_S2_25_source_crosswalk_20261118.json). The TTF frozen main is unchanged.
