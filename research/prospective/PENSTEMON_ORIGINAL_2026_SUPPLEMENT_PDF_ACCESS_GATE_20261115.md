# Penstemon 2026 phenotyped accession source: supplemental PDF access gate

**Registration time:** 2026-11-15. The 4-origin BioSample-only gate has already passed on GitHub Actions #43711234306, but none of the 34 measured flower specimens is yet linked to the verified 52 genome samples. This source-only check is for the *2026 original paper's* reported Table S1, not to infer phenotyped status from genome accession.

## Source authority

Stone et al. 2026 bioRxiv, DOI [10.64898/2026.07.21.739878](https://doi.org/10.64898/2026.07.21.739878), archival PMC [PMC13419788](https://pmc.ncbi.nlm.nih.gov/articles/PMC13419788/). The paper says Table S1 in Supplement 1 lists **collection sites and vouchers for 52 samples**, while **34 accessions** have measured floral traits; authors say analysis-ready trait tables and code will be made available on Zenodo following publication. Table S1 inclusion alone **does not** establish which 34 were phenotyped.

Only a finite list of official-host source candidate URLs may be checked. The script verifies genuine PDF signature and emits **file size / SHA256 / page count / which pages mention Table S1 or Table S2**. It does *not* export Table S1 content or identify individual vouchers yet; no DNA, phenotype, RI or ancestry values are decoded.

If official supplemental PDF is accessible, a **separate** original specimen-key extraction and manual source-table verification may be attempted, distinguishing source specimen identifiers (CAW, BWS, LW, etc.) from NCBI BioSample SAMN keys. A specimen label is not identical to a BioSample accession merely because species names agree.

If all direct retrieval attempts are blocked, report `SUPPLEMENT_PDF_SOURCE_UNAVAILABLE`; no guess at 34 specimen IDs, no claim of a four-origin measured trait panel, and no change to frozen TTF main science.
