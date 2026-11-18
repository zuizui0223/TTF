# Penstemon × TTF/CREST: original molecular material and four-origin admission gate (2026-11-15)

**Result status:** SOURCE CONFIRMED for 52 public whole-genome BioSamples; source-only PRELIMINARY AUTHOR TREE confirmed at 45 tips. **The 52 new WGS samples, 34 phenotyped accessions, and 239-species 43-locus tree are NOT joined.** No genetic or floral transfer metric has been computed. The frozen TTF-Q main results are not changed.

## Author publications and their distinct grains

| Original publication/source | Valid source unit and evidence | Admissibility |
|---|---|---|
| Stone et al. 2026 bioRxiv DOI `10.64898/2026.07.21.739878` | Four independently inferred hummingbird-pollination syndrome origins within *Habroanthus*; 52 unique sampled populations/accessions for whole-genome resequencing, 34 accessions phenotyped | Independent **origin** is the holdout unit; 52 samples are not 52 independent origins |
| ENA/NCBI BioProject `PRJNA1479764` | **52 distinct run IDs and 52 distinct BioSample accessions**, source-run metadata read from the ENA public API; species and submitter-alias fields were blank in the `read_run` source report | Confirms 52 sequencing source records, **not** a 52×34 specimen join or a 52-species panel |
| `lcwheeler/penstemon_novaseq` public source, GitHub last updated **2023-06-02** | Original author CDS and genomic-window trees each have **45 unique tip IDs**, with the **same 45 IDs**. A separate `target_trios/groups.tsv` maps **18 of these 45 tips** to 8 literal species labels. See [machine receipt](penstemon_2023_public_45_tips_and_18_groupkeys_receipt_20261115.json) and successful source workflow #43708795462 | Valuable historical source, **not guaranteed to be the final 2026 52-accession genome tree**. Do not auto-equate CAW sample labels to 2026 SRA accessions |
| Wolfe et al. 2026 Dryad DOI `10.5061/dryad.g79cnp66z` | Version 4 published 2026-11-02, original archives for ~239 species and 43 nuclear loci; RAxML/Astral tree archive digests pinned in [machine receipt](penstemon_239_nuclear_source_only_receipt_20261115.json) | Official individual tree downloads return **HTTP 401** for anonymous API. No tree tips or accession crosswalk established from this source |

## Four independent hummingbird-lineage contrasts

Original 2026 text identifies four evolutionary origins:

| Hummingbird-lineage group | Bee-lineage sister comparator | Old 2023 `target_trios/groups.tsv` taxon names represented? |
|---|---|---|
| *P. barbatus* | *P. virgatus* | Both mapped (five barbatus; two virgatus) |
| *P. cardinalis* | *P. glaber* | Cardinalis mapped (one); **glaber NOT mapped** in the old **18-row subset** |
| *P. eatonii* | *P. laevis* | Both mapped (five eatonii; two laevis) |
| *P. labrosus* | *P. deaveri* and/or *P. putus* | Labrosus and putus mapped (one each); deaveri not in the 18-row map |

**Strict source gate:** **3 of 4** published sister-contrast group pairs are explicitly represented on both sides in the older 18-row mapping; the fourth lacks a mapped *P. glaber* comparator. This is *not* evidence *P. glaber* is absent from the older 45-tip tree; it is only a taxonomic-metadata coverage limitation. A confirmed four-origin train/test analysis must use final source accessions and equivalent phenotypes, not old partial labels.

## Molecular hypothesis and contrast that remains scientifically open

The original 2026 preprint reports:
- Strong multi-trait phenotypic convergence across four hummingbird syndromes;
- Different anthocyanidin routes to red (mainly pelargonidin in *P. barbatus*, *P. eatonii*, *P. labrosus*, versus cyanidin in *P. cardinalis*);
- Extensive genome-wide incomplete lineage sorting; no detected adaptive introgression fueling all four hummingbird origins;
- No comprehensive mapped causal trait loci across all four origins. Therefore ILS versus shared ancestral variation versus new mutations is not solved for the developmental loci.

**TTF-P:** Will a trait-change direction calibrated on three independent hummingbird-syndrome origins predict the fourth independently heldout origin? The response is a quantitative floral trait such as nectar volume, pigment spectrum, tube dimensions, not a binary syndrome label derived from those traits.

**CREST/TTF-G:** After reconstructing history from outcome-independent, voucher-matched nuclear loci, does history add predictive power over the current floral phenotype and environment? A shared lineage topology alone is not evidence of a historical effect or cause.

**TTF-Q gate:** Require independently published accession×phenotype keys, study-specific tree/sample concordance, all four origins with bee sister comparators, measured response availability in every heldout origin and a source-defined correct taxonomic root. Hold out entire independent origin including all populations/accessions; *never* randomly split sibling or same-species individuals. Four origins produce only four evolutionary test units, making inferential power limited even when model scores are calculable.

**Current status:** `BIOPROJECT_AVAILABLE_52`, `OLD_NUCLEAR_TREE_AVAILABLE_45`, `OLDER_TAXON_MAPPING_ONLY_18`, `LATEST_WIDE_GENUS_TREE_AUTH_REQUIRED`, `PHENOTYPE_SOURCE_MATCH_UNVERIFIED`, `TTF_P_G_RI_NOT_EVALUABLE`.

## Reproducibility

- [ENA 52-run, 52-BioSample original-metadata receipt](penstemon_habroanthus_ena_52_accession_source_receipt_20261115.json), source run [#43708164658](https://github.com/zuizui0223/TTF/actions/runs/43708164658).
- [Older author 45-tip exact-match / 18-group metadata receipt](penstemon_2023_public_45_tips_and_18_groupkeys_receipt_20261115.json), source run [#43708795462](https://github.com/zuizui0223/TTF/actions/runs/43708795462).
- [Dryad wide nuclear phylogeny source inventory](penstemon_239_nuclear_source_only_receipt_20261115.json), source run [#43701619489](https://github.com/zuizui0223/TTF/actions/runs/43701619489).
- Public ENA `returnFields?result=sample` was diagnosed as **text/plain TSV** with first header `columnId, description, type`, source-only [run #43709185336](https://github.com/zuizui0223/TTF/actions/runs/43709185336). Original sample aliases still require a successful new sample-level report; no unsupported taxon inference.

All added code is separately reviewed in Draft PR #159 and is independent from the frozen TTF-Q main scientific claims.
