# Costus 2022–2026: public source lineage, matched taxonomy and TTF-Q stopping result

**Date:** 2026-11-15. **Status:** the source-only audit is complete for the 2022/2025 linked tree–scent dataset; `ROOTED_CLADE_NOT_EVALUABLE` and `UNROOTED_BLOCK_NOT_EVALUABLE` under two frozen prospective geometry families. **No floral scent response intensity, compound richness or molecular mechanism was fitted**, and no genetic/speciation-transfer effect is claimed.

## 1. Source availability

- Kay & Grossenbacher (2022): Dryad DOI [10.7291/D1C39G](https://doi.org/10.7291/D1C39G), version 4 / ID 190546. Official individual-file downloads for the tree and observed pollinator CSV returned HTTP **401**, even though file metadata are public.
- Darragh, Kay & Ramírez (2025): [OSF 2ap4k](https://osf.io/2ap4k/), **public original data accessible and SHA-pinned**. The file `Floral Scent Analyses/chronogram_v8_2.tre` is byte-identical to the 2022 Dryad tree (`sha256:24f1cbfaa82ce4d33aad1ba80787cc611a6730aafb38f8244954c4d0d50ae7a1`). This legitimately solves the tree-source access blocker, **not** the direct observed-pollinator CSV provenance.
- The 2025 original OSF also contains `sampleinfo_floral.csv`, `sampleinfo.csv`, `means_species.csv`, `phylo_info.csv`, and nuclear gene-family phylogenies; verified file details in [nested OSF index](ttf_crest_costus_osf_nested_manifest_20261115.py).
- The separate 2026 Costus pollinator-efficiency Dryad [10.5061/dryad.931zcrk0r](https://doi.org/10.5061/dryad.931zcrk0r) has source files but official original download returned HTTP 401; the primary quantified field system compares **C. kuntzei and C. wilsonii**, not 30 independent Costus taxa.

## 2. Actual source keys and structure (not phenotypic outcomes)

| Original dataset | Source unit | Exact rows |
|---|---|---:|
| 2025 `sampleinfo.csv` | scented floral plus source ambient samples | **168** (only 100 floral) |
| 2025 `sampleinfo_floral.csv` | floral samples | **100**, 30 *distinct original sample taxon codes* |
| 2025 `means_species.csv` | original phylogenetic specimen/group mean | **27** |
| 2025 `phylo_info.csv` | source phylogenetic `Name` crosswalk | **27** |
| 2022 `chronogram_v8_2.tre` | original tree tips | **69** |

The 27 original `means_species.csv` first-column labels match the 27 `phylo_info.csv` `Name` keys **exactly**, and all **27** occur among the 69 original chronogram tips (different from the 2025 gene/TPS subset with four nominal species). The original `phylo_info` source category is 14 Bee / 13 Hummingbird.

**Important:** 27 tips are **not independently verified 27 biological species**, because multiple source accessions or infraspecific forms occur (e.g., two separate `Costus_scaber_...` vouchers, *C. guanaiensis* varieties). The 100 floral samples use 30 shorthand names; numerical scent result values have not been opened. No automatic taxonomic synonym resolution. Pollination categories are inherited from source labels; some original labels were inferred from floral morphology and therefore **must not be called 27 directly observed pollinator guilds**.

## 3. TTF-Q source-only gates and terminal decisions

- **Rooted clade gate:** original tree parser reports `is_rooted = False`. Under preregistered rooted-only eligibility, [Actions #43700817708](https://github.com/zuizui0223/TTF/actions/runs/43700817708) returned `NOT_EVALUABLE_ORIGINAL_CHRONOGRAM_NOT_ROOTED`. The tree was **not** rerooted by assumption.
- **Separately documented, pre-outcome unrooted-neighborhood amendment:** [preregistered K4/K3/K2 design](COSTUS_UNROOTED_PHYLO_BLOCK_SOURCE_GATE_20261115.md), based on unrooted tip-to-tip patristic distances, average linkage and a fixed minimum of three test taxa with at least one source Bee and one source Hummingbird per heldout cluster. [Actions #43701139981](https://github.com/zuizui0223/TTF/actions/runs/43701139981) succeeded technically; *all three* K candidates were **not evaluable**.
  - K=4: groups of 18 (Bee12/Hummingbird6), 5 (Bee1/Hummingbird4), 3 (Hummingbird3 only), 1 (Bee1 only).
  - K=3: groups of 21 (Bee12/Hummingbird9), 5 (Bee1/Hummingbird4), 1 (Bee1 only).
  - K=2: groups of 26 (Bee13/Hummingbird13), 1 (Bee1 only).
- **Consequence:** neither an independently reconstructed pollinator shift nor a robust tree-blocked out-of-group scent prediction can be claimed under these designs. No K/subgroup threshold was tuned after seeing scent outcomes, and the original sources are **not evidence for an ecological null**.

Machine receipts:
[27 × 69 original tree match](costus_tree_2025_pollination_source_geometry_receipt_20261115.json),
[30/27 exact sample/mean name crosswalk](costus_osf_scent_27_means_30_sample_keys_receipt_20261115.json),
[unrooted source STOP](costus_unrooted_phylogenetic_neighborhood_stop_receipt_20261115.json).

## 4. Practical external system advancement

For genuinely replicated molecular evolutionary transfers, a **different source panel** is needed rather than a post-hoc rearrangement of Costus groups. The Penstemon *Habroanthus* 2026 preprint [doi:10.64898/2026.07.21.739878](https://doi.org/10.64898/2026.07.21.739878) independently reconstructs four hummingbird-syndrome origins with 52 genome-sequenced accessions but floral phenotypes on only 34 accessions; author claims at this stage are primarily morphological convergence and ILS, **not proof that adaptive genetic alleles are shared**. A newly published, **separate** Penstemon nuclear phylogeny covering 239 species and 43 nuclear loci ([Dryad 10.5061/dryad.g79cnp66z](https://doi.org/10.5061/dryad.g79cnp66z), Nov 2, 2026) offers a much broader prospective ancestry backbone. The 239-species phylogeny and 34-accession flower measurements must be crosswalked on strict original taxon and accession IDs before treating them as a common multi-level data panel.

**TTF/CREST claim boundary:** this report documents *where a historical phylogenetic state can be established from source* and *where the observational geometry fails*. It is not a test of an observed hidden-history effect, reproductive isolation, or universal transfer. Frozen prior TTF-Q, COI, He, Fst, breaks and P0/P1 results remain untouched.
