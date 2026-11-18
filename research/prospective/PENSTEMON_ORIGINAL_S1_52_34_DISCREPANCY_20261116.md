# Penstemon 2026 Table S1 mismatch: independently verified original-source HOLD (2026-11-16)

**This is a source-quality diagnosis, not a biological response analysis.** The original Europe PMC supplement `media-1.pdf`, SHA256 `3682720dad4f0621bebef3be7014aa667067c5f0dce9160ca8f4c4a758e7f943`, was checked on 2026-11-16 via [independent source-only Actions run #43874069838](https://github.com/zuizui0223/TTF/actions/runs/43874069838) using two independent PDF text extraction engines, including page 12 only **before** Table S2. Aggregate artifact ID `14228504100`; ZIP SHA256 `70ef924b7a741bdda1d968d01401488248aac36a94a48783da08cfbc6d3a9887`.

## What can be verified without inferring a missing source record

| Original-source evidence | Result | Admission |
|---|---:|---|
| Paper reports WGS source samples | 52 | Independently verified 52 NCBI BioSamples |
| Paper reports flower-phenotyped accessions | 34 | Not yet matched to voucher IDs |
| pypdf plain **complete** Table S1 voucher lines | **51** | Incomplete against the reported 52 |
| PyMuPDF sorted **complete** Table S1 voucher lines | **51** | Independent source-engine confirmation |
| Complete S1 lines marked `Phenotyped=yes` | **30** | Cannot equate to paper's 34 |
| Complete S1 lines marked `Phenotyped=no` | 21 | Source flag only, not proof of no partial field measurements |
| Voucher IDs not mapped to original NCBI taxon strings | `CAW-190` | **Taxon mapping HOLD**; do not invent synonym |
| NCBI 2026 taxon counts not equal to exact source PDF taxon matches | *P. cardinalis* NCBI 2 vs PDF 1; *P. comarrhenus* NCBI 2 vs PDF 1 | Missing/mismatched source label(s) unresolved |

**Why counting all ID occurrences is invalid:** the original PDF contains additional ID mentions and shorthand taxon labels, e.g., species-name prefixes paired with numbers. Across the original source page, loose token matches greatly exceed complete S1 voucher records. Counting all tokens as original specimens would falsely claim completeness. Sorted PyMuPDF extracts the **same 51 complete** voucher rows and the same 30 yes flags as pypdf; unsorted PyMuPDF separates PDF columns and cannot be parsed as rowwise records.

**Why the counts need not be exact biological contradiction:** The paper's four accessions measured in the field for a subset of flower traits could relate to the 30 vs 34 gap, but this is **only a possibility**; no correspondence is established by the source. It would be incorrect to reclassify four `no` records as `yes` without direct original specimen evidence. One missing *P. cardinalis* comparison and a `CAW-190` taxon string that differs from NCBI could reflect manuscript/metadata differences; **do not** decide equivalence from geography or morphology.

## Correct current TTF-Q disposition

- **PASS:** NCBI whole-genome study accession metadata 52/52; all four hummingbird-origin species and their required bee comparator species have genome records. No accessions/SAMN were inferred from row order.
- **HOLD:** published voucher × 34 phenotyped identity table, specimen-level link to SAMN, one unresolved original Table S1 voucher row and one unmatched taxon string.
- **HOLD:** true four-origin heldout trait prediction, CREST historical-genetic increment and reproductive isolation. No quantitative floral or genetic response was fitted.
- **NEXT source distinction:** independently index original *Table S3* sample identifiers for the published phenotype source, then visually inspect the authenticated Table S1 PDF page. If those original sources do not establish exactly 34 measurable accessions and strict voucher links, report `PHENOTYPE_SOURCE_NOT_MATCHABLE`; do not relax 52/34 admission to 51/30.
- The older 2023 45-tip tree is not directly equated to 2026 whole-genome BioSamples; the 239-species 2026 Dryad nuclear tree remains authenticated metadata with original tree archives HTTP 401.

## Provenance and files

- [52 BioSamples taxon source receipt](penstemon_habroanthus_52_ncbi_biosample_taxon_receipt_20261115.json)
- [Original 52-voucher prospective source rule](PENSTEMON_FULL_S1_PHENOTYPED_VOUCHER_SOURCE_RULE_20261115.md)
- [Original PDF two-parser identity-only diagnostic code](penstemon_original_S1_two_parser_52_51_source_diagnostic_20261116.py)
- [Original-source two-parser Actions run](https://github.com/zuizui0223/TTF/actions/runs/43874069838)

**Freeze:** Does not change submitted TTF-Q science or Senecio P0/P1 and does not assert a new biological transfer or null.
