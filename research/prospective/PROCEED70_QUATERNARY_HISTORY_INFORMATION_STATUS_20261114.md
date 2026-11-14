# PROCEED 70 plant species — Quaternary climate history information qualification

**As of 2026-11-14.** This document records **predictor geometry only**, not evolutionary rate, plasticity, adaptation, fitness, demographic genetic differentiation, or speciation.

## What actually ran

The source-only original-pair climate extraction [GitHub Actions #43221686424](https://github.com/zuizui0223/TTF/actions/runs/43221686424) recovered the exact SHA256-locked 21ka and 0BP CHELSA-TraCE21k BIO1/7/12/15 bands at **253 original source collection sites**, with **all eight values finite at every site**. Input rows: 1,561 (distinct from 1,580 pairs with both coordinates before original binomial/study eligibility), **70 plant taxa, 57 genera, 104 source studies**, completely study-and-species-disjoint train/evaluation panels **35/35 taxa**.

The precommitted response-blind source information test [GitHub Actions #43576786261](https://github.com/zuizui0223/TTF/actions/runs/43576786261) collapsed 1,561 original trait/effect rows to **187 unique study–species–ordered-geographic-site pairs** (102 train / 85 heldout). Every training species had equal total weight; within species, predictor source pairs carried equal subweight.

**Focal historical predictor:** for each source site, Euclidean magnitude of four local 21ka-minus-0BP climate differences standardized using training-species-equal current-climate SDs; the pair average of two site magnitudes is the single proposed historical exposure.

**Response-blind baseline:** four current-climate pair means, four absolute current-climate pair differences, geographic pair midpoint latitude/longitude, and great-circle distance (11 variables). In the 35-source-species training geometry, the weighted fraction of total historical-predictor variation not linearly reconstructed from this baseline was:

**0.5419597099656406 (54.20%).** The predeclared minimum source-only nonredundancy threshold of 0.05 was satisfied.

The predictor support was **not broadly uniform**: residual historical squared signal had species effective sample size **11.613 / 35**, and its largest training species contributed **20.923%**. The heldout training-projection squared error ratio was **0.659632**; this is NOT a heldout biological R², predictor-unique variance fraction, or evidence that history predicts evolutionary rates. The original 8-band numeric matrix SHA256 was reproduced exactly: `7d03b2f158ade223b752a2b7d6ac57c38fd163e0b44f5566a01e8e6fe6d6b4fb`. Original coordinate SHA256: `aed6486ee9aa49b0d454426c4d5964ebd79a8a96db2aac452c80fd5b7f52af32`.

## Scientific interpretation boundary

**SUPPORTED:** A late-Quaternary source-climate displacement magnitude contains distinct measurable variation beyond a specified *linear* contemporary climate/geographic baseline in a small, study-and-species-separated plant panel.

**NOT SUPPORTED YET:** That historical variation predicts genetically based quantitative trait change, reciprocal-transplant fitness, local adaptation, reaction-norm plasticity, reproductive isolation, or species diversification. The 70 species were eligible on geography and PROCEED's original `genphen=Genetic` study label; they were **not selected from demonstrated adaptation**. PROCEED's original publication explains that "Genetic" represents common-garden or quantitative-genetic efforts to account for heritable components, which is not a demonstration of adaptive benefit.

**LIMITATIONS:** (a) 187 independent predictor *geometries*, not 1,561 independent sites; (b) residual source signal concentrated in about 11.6 effective species; (c) present source coordinates need not reflect LGM population occupancy; (d) the 54.2% value is based on linear nuisance adjustment and does not exclude nonlinear climatic/geographic redundancy; (e) only 35 test species available, making broad claims sensitive to taxa and studies. No original trait-rate response, genetic variance, raw phenotype, Darwin/Haldane value, or fitness endpoint was opened or fitted in the source-only qualification.

## Frozen next gate

A second original-source replay under the same Ubuntu 24.04, Python 3.12.15 and exact pinned numerical/GIS versions has been initiated as [Actions #43577189439](https://github.com/zuizui0223/TTF/actions/runs/43577189439). It checks exact original coordinate and 8-column climate SHAs **and** equality of the 54.19597% residual information, residual effective species support and maximum contributor share. No biological-response columns are read.

Only after verifying this input repeatability should a new, independent and bounded **rate-of-heritable-phenotypic-change** endpoint be specified (including original PROCEED transformation/trait units and original study dependence). A positive result for a rate metric must still not be labeled adaptation without documented fitness/selection evidence. A negative result must not be tuned into a positive by re-opening old TTF-Q/COI/He/Fst/Jensen outcomes.

Source contracts: `proceed70_plant_history_information_contract_v01.json`; original site support `proceed_plant_chelsa_source_support_contract_v01.json`; original species/study source support `proceed_v62_plant_genetic_geo_contract_v01.json`. Source-only scalar receipt: `proceed70_history_information_first_pass_receipt_v01.json`. Draft PR #154; no change to `main`.
