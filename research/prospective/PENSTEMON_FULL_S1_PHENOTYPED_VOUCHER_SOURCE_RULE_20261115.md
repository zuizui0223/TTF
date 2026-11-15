# Penstemon 2026 Table S1 phenotyped specimen-key extraction rule

**Source-only registration before extracting all Table S1 rows (2026-11-15).** Original Supplement 1 PDF was verified twice from Europe PMC, 21 pages, PDF SHA256 `3682720dad4f0621bebef3be7014aa667067c5f0dce9160ca8f4c4a758e7f943`, Table S1 on PDF page 11. The limited first-20-line header preview (Actions #43712481678) explicitly reports a `Phenotyped` column with `yes`/`no` values alongside original population IDs (`CAW-`, `BWS-`, `LW-`, etc.).

This source-only operation is *not* opening floral phenotype magnitudes or SNP/fitness values. The `yes/no` indicator is missingness/availability provenance only.

## Predeclared conservative extraction

1. Fetch and check the same original Europe PMC PDF SHA256, then extract **only page 11 (1-based)** via the plain text extractor, which preserved the Table S1 title; layout mode demonstrably omitted it.
2. Accept exactly one source record per visible author voucher line beginning with an uppercase-letter prefix plus hyphen and 2–6 digits, optionally followed by an asterisk. The original Table S1 asterisk refers to the 5-taxon introgression analysis, **not phenotyped status**.
3. Require the line to end in `yes` or `no` for the source `Phenotyped` column; never impute an answer from the presence of an abbreviation in Figure 1.
4. Check exactly **52 distinct voucher IDs** and **34 `yes` statuses** against the original article's reported sequenced and phenotyped counts. Any mismatch is a **SOURCE_TABLE_FORMAT_OR_PROVENANCE_HOLD** and must not be silently repaired by dropping lines.
5. For provenance-only species labels, match an original PDF line's `P. speciesname` after removing intraword whitespace to exactly one of the **21 already authenticated literal BioSample taxon labels**. No synonym merging (e.g. `comarrhenus` vs `strictiformis`) or invented pair mapping. Report any unmatched or ambiguous rows.
6. Aggregate by four pre-fixed independent origins (derived × all bee sisters); note presence/absence of measured specimens separately from genome source. Retain source voucher IDs as published original metadata, not `SAMN` accessions, until an original-source specimen-voucher mapping is authenticated.

The code must NOT read Table S3 flower trait values, the PDF's quantitative physiology, or NCBI genotype/sequence. Success can establish that the four origins have **phenotyped specimens by published species labels**; it does not establish same-genome-individual causality or a TTF out-of-origin heldout prediction.

No changes to frozen TTF-Q main, original Senecio P0/P1 or Costus/Timema stops.
