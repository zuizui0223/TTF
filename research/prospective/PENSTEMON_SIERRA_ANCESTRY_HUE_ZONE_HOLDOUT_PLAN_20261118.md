# Independent Penstemon Sierra hybrid-zone × TTF source-locked exploratory pilot (2026-11-18)

**Design fixed before first full source-data phenotype/ancestry-value evaluation on this branch.** An earlier source preflight inspected only 4 individual rows and the complete metadata columns, and established counts: 421 distinct sample IDs, 5 population labels, 2 hybrid-zone groups and 357 rows containing ancestry coordinates. The pilot is accordingly **exploratory, not blinded confirmatory**. It is distinct from PR #159, which evaluates the four independent *Habroanthus* hummingbird-syndrome origins and currently lacks a 34-phenotyped-voucher mapping.

## Original independently published empirical source

Stone et al. (2026), *Evolution*, online 2026-11-05, **DOI 10.1093/evolut/qpag199** (superseding the 2026 bioRxiv preprint 10.64898/2026.01.17.700129). Public author repository [benstemon/sierra-hybrid-zone](https://github.com/benstemon/sierra-hybrid-zone), original source commit `0ba4842824bcd9ddc824eedea303ad1c176fc82b`, table `1.phenotype_curation/data/allcombined_phenotype_admixture_subpopulations_spreadsheet-sierras.csv`, original Git blob `7cd5abd0f2a5ad6de1563d3ddd2c0d857b037de8`.

**Original grain:** one individual/sample `newID`; `modified_hue` is the original author's continuous flower hue outcome, `K1` a source-published genome-wide ancestry proportion estimated from the original *combined* genotyped sample panel. `subpopulation` distinguishes `hybrids_virginialakes` and `hybrids_gemlakes`, each containing one set of the same *P. davidsonii* × *P. newberryi* hybrid-zone samples. A source-level count showed 63 / 220 hybrid rows in these two zones (but not whether every one has K1 and hue). No other `region`, parent taxon, pooled GWAS-significant marker, or gene-based classifier is used.

## Frozen single predictive question

Given **already estimated global ancestry K1**, is the linear *ancestry-to-floral-hue relationship* calibrated on hybrids from one of two independently sampled Sierra hybrid zones predictive for the entirely withheld other zone?

This is a **same-species-pair, cross-zone pilot**, NOT a cross-species TTF result, not four independent speciation events, and not evidence that neutral history adds information beyond current genotype. It is an operational predictive test of one already available molecular-to-phenotype interface.

### Exact response and population eligibility

- Only rows with literal `class == hybrid` and `subpopulation` in `hybrids_virginialakes` / `hybrids_gemlakes`. Parent samples and 64 unassigned/NA rows are excluded **without inspecting hue values**.
- Retain individuals with a finite numeric `modified_hue` and `K1` (range K1 ∈ [0,1]). Require a unique `newID` overall and within zones. Do not impute phenotypes, ancestry or zone; all missingness is reported per zone.
- Train-only ordinary least squares `modified_hue = a + b*K1`, with intercept. If training K1 variance is zero, return NOT_EVALUABLE. No transformation beyond the original author's `modified_hue`; no regression regularization or hyperparameter search.
- Exactly two folds: train `hybrids_virginialakes`, evaluate `hybrids_gemlakes`; then train `hybrids_gemlakes`, evaluate `hybrids_virginialakes`. Require ≥20 eligible hybrid individuals in each zone; do not use heldout `modified_hue` to select/flip gene direction or tune model.
- Fixed benchmark is **training-zone mean hue** for every heldout individual. Compute each fold's holdout RMSE, MAE, squared-error skill `1−SSE(OLS)/SSE(train-mean benchmark)`, and heldout Pearson correlation (descriptive, not an independent fitting objective). Report source observation counts and both fitted training coefficients; do not pool individual rows across the two folds to manufacture independent ecological replication.
- Report the two fold results separately and the pair-equal mean skill **descriptively only**. Two zones provide no population-level predictive interval or hypothesis-test p-value; no claim of generality to another species pair.
- The original K1 was inferred jointly with test-zone genotypes, hence **this is transductive ancestry-feature transfer**, not strictly fully independent prospective genotype prediction. Also phenotypes are observational; ancestry–hue association is not proof of a causal color locus or reproductive barrier.
- If a failure gate triggers, stop without re-tuning, and report NOT_EVALUABLE rather than silently changing the population subset.

### Publication novelty ceiling

The published article already reports genome/flower divergence can decouple from genomic reproductive barriers, including a conspicuous hue-associated region with relatively permeable clines. We do **not** claim discovery of that mechanism. The exploratory methodological addition is only the *source-pinned out-of-zone fixed-coefficient ancestry→hue predictive comparison* and its interpretation under TTF's transfer-versus-history distinctions.

No frozen TTF-Q, TTF initial genetic outcomes, Penstemon PR #159, or original Senecio results are changed. This branch requires a separate draft PR for review; no automatic main merge.
