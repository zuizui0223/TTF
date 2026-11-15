# Can Quaternary climate history predict speciation, adaptation, or plasticity?

**Source-first feasibility and prospectively bounded programme. 12 November 2026.**

This is **a separate research programme** on branch `prospective/quaternary-phenotype-source-screen-20261112`, based on the **published metadata of prospective external responses**. It is not an amendment of the closed TTF-Q, post-IBD COI, He, Fst or Jensen genetic-break experiments. No new association coefficient, p-value or predictive gain has been calculated.

## Answer and current decision

**Yes in principle, but the three phenomena demand different observable responses and time scales.** TTF-Q establishes that a late-Quaternary climate relation is not fully redundant with modern climate and geographic opportunity; it does NOT imply that the extra information predicts any evolutionary outcome.

| Endpoint | Required biological response | Public source currently identified | Current source-only decision |
|---|---|---|---|
| **Phenotypic plasticity** (near-term candidate) | A within-genotype/population **reaction norm**, e.g. acclimation response ARR = (thermal limit at treatment 2 − treatment 1)/(treatment temperature 2 − treatment 1); CTmax and CTmin treated separately | Weaving et al. 2022 OSF 102 insect species and 1,374 effect sizes; Barley et al. 2021 source-labelled 19 animal species, 78 populations | Public effect tables exist. Weaving's verified source XLSX has standardized latitude and latitude bin but **no raw source longitude**; its exact local LGM climate exposure is not identifiable from that spreadsheet alone. Barley has only 12 terrestrially georeferenced animal species passing descriptive ≥2-site census; six source-population coordinate keys conflict. **NOT EVALUABLE as a macro test today.** |
| **Heritable/local adaptation** (next metadata candidate) | Common-garden among-genotype/population differences or reciprocal transplant **fitness-by-origin × site** (home–away advantage); wild phenotype change alone cannot establish adaptation | PROCEED v6.1: 9,263 records, 326 original studies, 428 taxa; includes some field phenotypic and some common-garden/quantitative-genetic designs | Identify original-study heritable-design flag, original coordinates and population overlap **before** checking trait effects. Dataset size does not equal a clean genetically based 428-species panel. |
| **Speciation** (different macroevolutionary programme) | Dated sister-lineage splits with evidence of reproductive isolation/lineage persistence, not simply current genetic distance | VertLife mammal phylogenetic posterior ~5,911 species; species-level phylogeny with some imputed taxon placements | A 21-ka climate snapshot **cannot be assumed to cause branches dated hundreds of thousands or millions of years earlier**. Need longer palaeoclimate time windows and ancestral-range history; phylogenetic diversification rates alone are not proof of speciation mechanism. |

## Verified source-only results

1. **Barley**: original `jmbarley1/plasticity` Git commit `c1416c8925d6f1509b94fb61a9e7973f4cf03cb4`, raw CSV blob SHA1 `8e960f1e43f8d4f4ca4c9521c41d8576e97beaec`. A pre-phenotype census found 326 source rows, of which 210 were tagged terrestrial; 22 terrestrial rows represented non-animal *Arabidopsis* records, leaving 188 terrestrial-animal rows (Arthropoda 149; Chordata 39). Twelve animal species and 12 original studies had two nominal source populations; 10 genera. **Six conflicting origin-coordinate keys, all assigned to *Temnothorax curvispinosus*, remain unresolved**. Species counts and the publication's broader 19-species sample cannot be assumed identical. Nothing is inferred about ARR direction or LGM exposure.
2. **Weaving**: public OSF node `cbhv4`; verified exact `Supplementary Data 1.xlsx` 351,661 bytes SHA256 `73b2d9be32d54177e009c2f434dcebd1d918282a12d9a52aaf3a970c96f6ff4e`. Main data worksheet `A1:AR1377` (declared Excel dimension, **not a verified row count**), header includes `studyid`, `popid`, `speciesid`, `full_sp`, `stand_lat`, `lat_bin`, `ARR`, `ctlimit` and design variables. **No exact unstandardized latitude/longitude pair in source headers.** An explicit original-study location reconciliation is required for a population-exposure estimand.
3. **PROCEED**: its authors separate phenotypic-only field changes from common-garden or quantitative-genetic studies. Even a genuinely heritable trait difference does not automatically establish *adaptive fitness benefit*. Source-only metadata discovery is the active next gate; detailed per-record heritable/field classification and locality support are not yet authorized.

## Candidate ecological mechanisms to predeclare on a NEW response-disjoint design

### A. Plasticity versus historical environmental predictability

Predictors should distinguish: (i) late-Quaternary environmental displacement magnitude; (ii) persistence of suitable local climate/thermal refugia; (iii) current and palaeo spatial niche heterogeneity; and (iv) **environmental reliability at the organism's generation/acclimation timescale** (requires different meteorological data, not automatically inferable from 21ka/0BP rasters). Select **one** main historical contrast and a bounded set of secondary mechanistic predictions *before* observing the new ARR association. Climate history can select for canalization, bet hedging or plasticity depending on cue reliability, costs, and temporal predictability, so no universal positive slope can be assumed.

### B. Genetic adaptation versus ecological response

If one can isolate common-garden or reciprocal-transplant studies, test whether historical refuge persistence/temporal climate contrast predicts **heritable differences in trait optima**, or home-versus-away fitness interaction. Control contemporary climate, habitat, dispersal, phylogeny, study and assay traits. Treat genetic adaptation as **not evaluated** if only wild field phenotypes or standardized trait changes are recorded.

### C. Speciation with compatible chronological scope

Map the dated late-Pleistocene and Holocene window only to divergences with a credible interval plausibly spanning that period. A broad lineage-diversification test must use longer Pleistocene/Pliocene environmental histories and the posterior uncertainty of tree topology and branch dates. If contemporary distributions are used as glacial ancestral occupancy surrogates, report that as a separate, noncausal comparison.

## Source-blind macro qualification and abstention rules

- Require **≥50 terrestrial independently identified taxa**, ≥15 source studies, ≥10 original genera, and ≥20 taxa per disjoint training and evaluation arm for the first exploratory held-out plasticity/adaptation benchmark; this is a self-imposed *design minimum*, not evidence of adequate power. Recalculate design-specific power after source eligibility, prior to any association opening.
- Exactly georeferenced **original experimental collection populations** are the primary exposure units. Species-wide GBIF/SDM LGM environmental-envelope measures are alternative **species-range proxies**, not replacements silently attached to the same response estimand.
- Every prospective taxon and study must belong entirely to its designated evaluation arm. Related effects within a publication/population/lineage must be clustered or blocked. Bootstrap units are source studies/species, not 1,374 independent effect sizes.
- Controls: modern current climate, seasonality, altitude/latitude, macroregion, source-study methods, acclimation temperature contrast and duration, premeasured basal thermal tolerance, taxonomic phylogeny and sampling intensity; report overlap, collinearity and spatial leakage.
- **Predictive criterion**: compare declared source-only controls against one a priori history increment on withheld studies and species, with interval uncertainty and a relevant control for observational geography. Distinguish **predictor uniqueness** from biological held-out gain.
- Existing TTF results remain frozen; source census cannot be called negative adaptive/plasticity evidence. If the chosen external panel fails source support, mark `NOT_EVALUABLE` instead of changing endpoints, inferring missing origin sites from GBIF centroids, or tuning to a positive.
- Cross-taxon conclusions from a single 12-species or highly Diptera-skewed panel require explicit limitations. A larger compiled dataset is valuable only if biological endpoint and geographic provenance genuinely qualify.

## Immediate routing

**Priority 1:** Complete source-only PROCEED v6.1 file catalog and identify whether original per-record heritable-design flag and collection site coordinates are available without reading trait change values.

**Priority 2:** Recover experimental origin locations for the Weaving species/study/population IDs *as an independent source document with exact provenance*; if most cannot be recovered, either use a separately declared **species-range historical niche** test or do not make a LGM–plasticity prediction.

**Later:** Time-calibrated speciation/late-Pleistocene divergence using a genuinely independent dated tree and predeclared ancestral spatial history; do not use the frozen genetic negative outcomes as a new tuning target.

Source identity receipts: `barley_2021_source_only_pilot_receipt_v02.json`, `weaving_2022_source_header_receipt_v01.json`, and `proceed_v61_source_metadata_contract_v01.json`. No experimental response values were inspected in these feasibility checks.
