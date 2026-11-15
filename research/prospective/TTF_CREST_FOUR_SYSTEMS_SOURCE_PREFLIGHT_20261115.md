# CREST × TTF: independent evolutionary systems, source-only preflight (2026-11-15)

**Status: SOURCE-ONLY / EXPLORATORY. No biological response, genotype effect, fitness effect, or transfer score has been newly fitted or opened.** This is a separate candidate-screening branch; it does not change the frozen TTF-Q, COI, He, Fst, or break-location claims on main.

## Primary objective

Can a relation learned in one lineage, habitat pair, island, or pollination transition predict a truly withheld unit after correcting for morphology, geography, historical ancestry, and sampling? Compare response grains separately: morphological trajectory, genomic mechanism, and reproductive isolation. A successful internal analysis is not cross-species transfer unless the unit held out is a species.

### Source audit matrix

| System | Directly verified source and grain | Possible TTF held-out unit | Current stopping issue |
|---|---|---|---|
| **Senecio lautus** | Author GitHub at pinned commit below: 605 phenotype individuals in 20 populations; 481 individuals in 16 populations, forming **8 distinct habitat pairs** in paired file; pair-level means in 16 rows. A separate genotype analysis has **9** Dune–Headland pairs and 9,269 SNPs. | Locality pair (TTF-P); later independent lineages only if cross-species support is added | **Individual-level genotype/phenotype join NOT verified.** Only eight of nine named genotype pairs occur in paired morphology. A genome outlier list made using the full data is not an independent test predictor on held-out pairs. |
| **Timema** | Dryad 10.5061/dryad.q573n5tpk describes **42 host-isolation pair records** (8,749 preference trials) and **30 pair-level mating datasets** (2,436 trials), matched by field `Pair`. Includes divergence-time TSV, nuclear alignment, and six species-specific variant files. Paper covers nine sexual species. | Species-disjoint, population-disjoint host/mating pair evaluation (TTF-RI) | Raw file headers/keys and graph components **not yet locally fetched or audited**. Pair endpoints recur; never treat all 42 dyads as independent units. |
| **Costus** | Dryad 10.7291/D1C39G has `chronogram_v8_2.tre`, `Pollinator_observations.csv`, and species morphology files; Darragh et al. 2025 assayed **100 scent samples / 30 species**, including eight reported independent shifts; analysis scripts and scent data are deposited on OSF 2ap4k; RNA-seq PRJNA1010186 examines a **restricted three-species clade**, not all 30 species; assembly PRJNA1091680. | Independent pollinator-transition/clade holdout (TTF-P) | 2022 tree tip / pollination label / 2025 scent sample matching **not audited**, and pollinator labels predicted from morphological traits must not be used as independent predictor truth. Field versus greenhouse origin is a potential source confound. |
| **Metrosideros** | Choi et al. 2021: **131 Hawaiian individuals / 11 taxa**, 22,511,205 variant sites; NCBI PRJNA534153 and Zenodo 4264399. Ten taxa are sampled as one population each; the widespread glaberrima taxon has three populations. | An island- or lineage-held-out genomic model only after comparable trait/site support | Taxon and locality are heavily confounded; no matched independent phenotypic-response table was verified. **HOLD** for a confirmatory TTF test. |

## Direct author-data audit: Senecio

Source repository: https://github.com/OB-lab/James_et_al._2021-Evolution
Pinned commit: `33222710aca75f79153e4118b409e7d77787d045`

| Original path | Blob SHA | Verified metadata |
|---|---|---|
| `phenotype/data_files/allTraits.txt` | `30f3edbdb486733306529e58135e6f1c205c0184` | 605 observations; 20 populations |
| `phenotype/data_files/allTraitsPairs.txt` | `ac8183378ea0778fa1a2fd42fd6a6d94d7008a85` | 481 observations; 16 populations; 8 named pairs |
| `phenotype/data_files/traitMeansPairs.txt` | `400f73bf21fefb6fcb775396f256c3447c755595` | 16 population means; 8 named pairs |
| `genotype/vcf_file/ESC_rel_50pp_80md_HWE_pairs_MAF0.05.vcf.tgz` | `fe1a452b183e0690e013cef57e59bedb13c7e072` | Binary VCF archive identified by author Git tree; genotype rows **not inspected** |

Morphology named pairs: `D00-H00`, `D01-H01`, `D02-H04`, `D03-H02`, `D04-H05`, `D05-H06`, `D12-H14`, `D32-H12`. Four populations not in paired morphology: `D09`, `D35`, `H03`, `H07`. Genotype analysis has an additional `D14-H15` pair, confirmed by a distinct author-file name; the other eight have same-named comparison files. This is a **population-pair name overlap**, not proof that individuals are shared or that all analysis outputs have common samples.

*The README reports population-mean imputation for approximately 11% of morphology samples and standardized phenotypic variables. These are transformed phenotypes, not independent raw field measurement columns.*

**Risk:** using published candidate SNP/outlier sets selected over all nine genomic pairs as train features in a leave-one-pair-out evaluation would leak test data. Proper test requires recomputing feature selection in each training fold from an admissible genotype matrix, or restricting the test to a precommitted phenotype-only ecological predictor.

**Unit restriction:** eight locality-pair comparisons, not 481 independent geographic experiments. A locality holdout is not a species holdout. A phylogenetic-clade holdout still needs a documented clade assignment and sample support.

Reproducible metadata-only check: `research/prospective/ttf_crest_four_systems_senecio_preflight_20261115.py` (fetches exact pinned bytes, checks Git blob SHAs, reads only population/paired-ecotype fields; fails on drift).

## Concrete next admission tests (not outcome analyses)

1. **Senecio:** confirm the eight pair IDs, geographic distinctness, matching genotype pair sample sheet and whether any independent ancestry tree permits grouping into at least two withheld clades. Do not reuse full-data SNP outlier selection on held-out locality.
2. **Timema:** obtain versioned Dryad `HostRIDat.csv` and `CountData.csv` with source checksum; reconcile `Pair` keys, explicit species/population endpoints, 42 versus 30 pair coverage, and connected components under shared biological endpoints. Determine whether truly disjoint species evaluation is possible **before** outcome values are decoded.
3. **Costus:** build strict source-only 2022 tree / measured pollinator / 2025 scent / OSF/SRA crosswalk. Derive independent shift/clade eligibility from tree without using scent as a selection criterion. RNA-seq is limited to a subset of species.
4. **Metrosideros:** obtain Zenodo file inventory and sample-to-taxon-to-island metadata; identify independent measured traits and multiple populations per lineage. Do not report TTF transfer when taxon predicts location by construction.

## Source identities and scope

- Senecio author data: https://github.com/OB-lab/James_et_al._2021-Evolution ; primary paper https://doi.org/10.1111/evo.14387 .
- Timema species-pair archive: https://doi.org/10.5061/dryad.q573n5tpk ; primary paper https://doi.org/10.1038/s41467-024-46294-9 ; SRA PRJNA356405.
- Costus phylogeny/morphology: https://doi.org/10.7291/D1C39G ; 2025 scent paper https://doi.org/10.1093/molbev/msaf027 ; OSF https://osf.io/2ap4k/ ; BioProjects PRJNA1010186 and PRJNA1091680.
- Metrosideros genome study https://doi.org/10.1073/pnas.2023801118 ; Zenodo https://doi.org/10.5281/zenodo.4264399 ; NCBI PRJNA534153.

**Evidence level:** Senecio author-source records were read directly with only key and header fields used. For the other three systems, source publisher/repository descriptions and paper methods were verified, but raw per-row records and cross-dataset joins were not inspected. No empirical transferability result is claimed. Biological novelty would require calibrated, disjoint held-out comparisons, not descriptive data inventory.
