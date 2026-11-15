# Penstemon × TTF/CREST: new independent empirical system (2026-11-15)

**Status:** source-only intake. No new molecular, floral-trait, adaptive or reproductive-isolation response has been opened or modeled. All frozen TTF main results are untouched.

## Why this independent system

1. Wolfe, Blischak, Stone & Kubatko (2026), *Phylogenetics of a rapid, continental radiation: Diversification, biogeography, and circumscription of the beardtongues (Penstemon; Plantaginaceae)*, Dryad **10.5061/dryad.g79cnp66z**, first released **2026-11-02**: **239 species / 43 nuclear loci**, ~84% genus coverage. Original archives listed as `1.RawFASTAs.zip`, `2.RAxML.zip`, `3.Astral.zip`, `4.DivergenceDating.zip` and related scripts / provenance.
2. Stone, Wheeler, Lambert, Williams & Wessinger (2026), *Rapid floral syndrome convergence in Penstemon through independent genetic variation*, bioRxiv DOI **10.64898/2026.07.21.739878**, published 2026-07-24: **four independent hummingbird syndrome origins** in section Habroanthus, **52 sequenced accessions representing 20 species and one outgroup**, and **34 floral phenotyped accessions**. Source reported major ILS, not significant adaptive introgression as repeated shift source. As of the preprint version reviewed 2026-11-15, code/data availability states release on Zenodo *following publication*, so an authenticated phenotype×nuclear voucher join is **not yet proven**.
3. Wessinger et al. (2023), PLOS Biology DOI 10.1371/journal.pbio.3002294: 229 individual genomes from a different Penstemon species complex with few (~20) species-diagnostic genomic regions. **Not** a 229-species dataset and not independent four-origin replication.
4. The separate 2026 Penstemon hybrid-zone paper in two species (P. davidsonii/P. newberryi) describes a large anthocyanin floral-hue locus with shallow genomic cline and other steep clines; this is a *possible future* barrier-vs-phenotype mechanistic external case, not proven to overlap with the 239-tip historical tree/vouchers.

## Immediate source-only task

- Determine the **exact** Dryad version and file IDs/digests for the Nov 2 nuclear-phylogeny authority. Fetch only metadata first.
- Test if the source RAxML/ASTRAL tree ZIPs are accessible without authentication; avoid a 231 MB all-data download, and never export genome reads.
- Separately document the 52 sequencing / 34 phenotype sampled figures and the 20 sampled in-group species plus one outgroup as different observation grains.
- Seek exact original source taxon/voucher keys. **Only after** a code/tree/trait strict join, check whether 4 derived pollination origins retain ≥2 different source-background lineages in both species-disjoint train and evaluation panels. No marker or subgroup may be selected based on the biological responses.

## Admission / claim firewall

- Primary user-facing biological question: whether the genetic mechanisms underlying similar derived floral syndromes predict flower traits and/or reproductive isolation in unseen evolutionary origins, *not* whether a guild label correctly classifies a flower in the same source.
- TTF-P (morphological prediction), TTF-G (genetic differentiation) and TTF-RI (isolation) are separate outcome definitions.
- A high phylogenetic tip count by itself does **not** prove matched floral traits or source-independent biological outcomes. Do not treat source-inferred pollination syndrome as field-measured pollinator visitation. Do not conflate 229 individual genomes with 239 species.
- If the trait original data are unavailable, report `PHENOTYPE_SOURCE_NOT_MATCHABLE`; if the source tree is unavailable, report `NUCLEAR_PHYLOGENY_RAW_UNAVAILABLE`. Both are **NOT_EVALUABLE**, not biological nulls.
- No changes to closed TTF-Q/COI/He/Fst/break-location endpoints or frozen results.
