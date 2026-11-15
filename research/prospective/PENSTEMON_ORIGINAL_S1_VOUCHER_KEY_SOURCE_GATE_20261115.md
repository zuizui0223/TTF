# Penstemon 2026 authentic Table S1 specimen ID inspection (2026-11-15)

Source-only escalation approved by the preceding source-gate results, before inspecting the actual supplementary Table S1 rows.

- Original article: Stone et al. (2026), bioRxiv doi:10.64898/2026.07.21.739878, PMC13419788.
- Authentic official Europe PMC `supplementaryFiles` ZIP, 9 members. The ZIP bytes are dynamic, so ZIP-container SHA is **not** considered immutable. Its inner `media-1.pdf` content is confirmed identical across two independently fetched official copies with SHA256 `3682720dad4f0621bebef3be7014aa667067c5f0dce9160ca8f4c4a758e7f943` and size 2,808,017 bytes.
- Original PDF has 21 pages. Page **11** (1-based) contains Table S1. Other headings: S2 page 12, S3 page 13, S4 page 14.

**Now inspect only page 11** to determine whether Table S1 exposes a reliable key connecting NCBI BioSample IDs and original author sample/voucher IDs and/or a field indicating which of 52 specimens have phenotypic measurements (34 reported in paper). Preserve exact original author ID strings. Print table title + column schema and a small source-layout excerpt *only*, never any floral numeric outcome, SNP genotype, sample DNA sequence, or reproductive isolation.

If Table S1 is a 52-specimen voucher inventory but lacks explicit measured phenotype membership or SAMN IDs, classify that join as `UNVERIFIED` and do not infer it from row order, common species name, or picture annotations. A Figure 2 point labeled with an accession does not automatically have a matching BioSample ID without original source-key evidence.

Next source-level inspection, if necessary, is a **separate** prospective read of exact original BioSample public identity fields (e.g. specimen_voucher / sample_name), not genetic or clinical sensitive metadata. Keep original 2023 tree CAW vouchers independent of new 2026 SAMN until verified.

No frozen TTF main scientific endpoint reopens and no biological effect-size estimate is approved.
