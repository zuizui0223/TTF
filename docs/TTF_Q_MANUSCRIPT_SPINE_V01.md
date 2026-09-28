# TTF-Q manuscript spine v0.1

Status: response-blind methods result assembled; relational v0.3 remains closed.

## Working title

**Evaluable envelopes for dyadic ecological predictors: response-blind decomposition of information, calibration and detectability**

Alternative ecology-facing title:

**Current and historical climate relations contain distinct pairwise ecological information beyond endpoint identity and geography**

## Core problem

Source-target ecological predictors are often treated as measured covariates and
tested directly against an outcome. In dyadic designs this can fail for at least
four different reasons that a single significance test cannot distinguish:

1. the relation itself is poorly measured;
2. endpoint identity or nuisance structure absorbs most relation variation;
3. usable information is concentrated in a few source or target taxa;
4. the intended estimator has adequate power only over part of the plausible
   effect-size / heterogeneity domain.

TTF-Q characterizes these properties before any biological response is opened.

## Methodological contribution

TTF-Q defines an **evaluable envelope** with distinct axes:

- exact source/target fixed-effect information survival;
- incremental information survival after declared controls;
- inverse-Herfindahl dyad/source/target signal concentration;
- null calibration under the intended dependence model;
- conditional effect-size × private-heterogeneity detectability;
- optional relation repeatability from response-blind resampling.

The key output is not PASS/FAIL. For each nuisance regime, TTF-Q first asks whether inference is calibrated. Only calibration-qualified regimes receive an evaluable MDE. This distinguishes non-identifiability, concentrated signal, anti-conservative inference and insufficient power rather than collapsing them into one stopping label.

## Empirical case-study ladder

All results below use external ecological data only. No nucleotide identity,
pairwise genetic distance, T_st, beta_R or beta_hist was opened.

### B1 — present-climate relation

On the shared 15,625-dyad development geometry:

- current-climate similarity retains 0.416 of its standardized variance after
  exact source/target identity;
- before geography, 0.406 of its total variance remains after endpoint identity,
  lineage and sampling-imbalance controls;
- after adding continuous 500-km geographic co-opportunity, 0.285 remains;
- geography therefore accounts for about 0.298 of the previously unique
  present-climate information;
- remaining signal is distributed across effective 61.9 sources and 74.6
  targets, with no source carrying more than 3.5% of squared residual signal.

This is a positive ecological-information result, not evidence for a genetic
transferability effect.

### C1 — historical climate beyond present climate

Across 653 climate-admissible species:

- raw historical/current relation correlation is moderate
  (Spearman 0.621 development; 0.592 confirmatory);
- after exact endpoint identity and controls including current climate,
  history-specific total unique variance is 0.365 development and 0.381
  confirmatory;
- effective signal support is broad: 78.9/81.8 source/target species in
  development and 139.1/151.9 in confirmatory.

Thus historical displacement similarity is not reducible to present climate,
lineage, sampling imbalance, or endpoint identity.

### C2 — geography decomposition

Adding continuous response-blind geographic co-opportunity reduces the
history-specific unique fraction:

- development: 0.365 -> 0.299; retained fraction 0.820;
- confirmatory: 0.381 -> 0.316; retained fraction 0.830.

Only about 17-18% of the C1 unique historical information is attributable to
this geographic co-opportunity measure. Roughly 82-83% survives.

## Empirical detectability result

Using the same two-way source/target dependence-aware estimator and the frozen
synthetic nuisance structure, TTF-Q reports the smallest effect-size grid point
whose Wilson 95% lower power bound reaches 0.80.

| private heterogeneity A | current climate | history given current |
|---:|---:|---:|
| 0 | 0.02 | 0.02 |
| 1 | 0.02 | 0.03 |
| 2 | 0.05 | 0.05 |
| 3 | 0.05 | 0.08 |

Null rejection rates remain near the intended one-sided alpha=0.025:

- current: 0.032, 0.027, 0.032, 0.018;
- historical: 0.022, 0.031, 0.023, 0.033.

The old single-cell reference of beta=0.03 at A=3 falls below the empirical
grid MDE for both relations on this shared geometry. That does **not** rewrite
the old relational-v0.3 decision. It demonstrates why a single worst-case
power gate is an impoverished summary of an otherwise informative design.

## Known-truth method validation

The decomposition is not justified only by the climate case study. A pre-result-frozen benchmark constructed dyadic predictors with exact truth for endpoint-retained fraction, baseline-control survival, and geographic survival. Across **243** scenarios, TTF-Q recovered all nested information fractions with maximum absolute error **3.33e-16**. Deliberate source and target concentration was detected in **100%** of the corresponding comparisons.

A second frozen benchmark compared four designs under the same legacy-style single cell, `beta=.03, A=3`. All four receive the same negative binary classification, but the calibration-qualified envelope separates the mechanisms. The broad reference design is calibrated with MDEs `.05/.08/.08/.10` across `A=0/1/2/3`. Endpoint loss and control redundancy remain calibrated but require larger effects (`.08/.10/.15/.20`). A source-concentrated design retains the same unique relation variance as the reference but has null rejection around **6–7%** at every tested `A`; it therefore receives **no evaluable MDE** under this estimator.

This benchmark establishes two points that a binary gate cannot represent: identifiable-information loss degrades detectability, while concentrated relation signal can instead invalidate calibration. Effective endpoint counts are therefore concentration diagnostics, not effective sample sizes.

## Main claims

1. **Dyadic ecological predictor quality is multidimensional.** Measurement,
   identifiability, signal concentration, null calibration and detectability are
   distinct properties and should not be collapsed into a single gate.
2. **An evaluable envelope is more informative than PASS/FAIL.** Known-truth
   designs with the same binary label can have different calibrated MDE
   profiles, while another design can fail because inference is anti-conservative
   rather than because information or power is low.
3. **Present and historical climate carry distinct relational information.**
   Present climate remains informative beyond endpoint/lineage/geography
   structure, and late-Quaternary displacement adds a second layer beyond
   present climate.
4. **The historical layer is not mainly geography in disguise.** About
   82-83% of C1 unique historical information survives explicit geographic
   co-opportunity adjustment in two independently assigned panels.

## Claims explicitly not made

- no genetic-transferability effect is estimated here;
- no v0.3 B or C confirmatory test is reopened;
- no old alpha is recycled;
- no NOT_EVALUABLE state is converted to a biological null or positive;
- synthetic detectability is not evidence that a biological effect exists;
- current and historical unique-variance fractions are predictor-information
  decompositions, not partitions of biological-response variance.

## Figure plan

**Figure 1 — TTF-Q architecture.**
Raw relation -> source/target FE -> nuisance controls -> geographic opportunity
-> information concentration -> detectability surface.

**Figure 2 — ecological information ladder.**
B1 current climate and C1/C2 historical relation, showing retained unique
fraction at each conditioning layer for development and confirmatory panels.

**Figure 3 — known-truth validation.**
Prescribed vs recovered information-survival fractions plus broad/concentrated
endpoint-support diagnostics.

**Figure 4 — calibrated detectability envelopes.**
Effect size on x, private heterogeneity A on y, with calibration-qualified
regions distinguished from anti-conservative regions; include the known-truth
single-gate comparison and empirical climate cases.

**Figure 5 — endpoint support.**
Source and target signal-share distributions, explicitly presented as
concentration diagnostics rather than effective sample sizes.

## Remaining optional extension

Response-blind occurrence-resampling repeatability can be added as a
measurement-error layer. It is not required for the now validated
information/calibration/detectability framework and should not delay the core
methods manuscript.

## Provenance anchors

- `benchmarks/frozen/ttf_q_c1_result_receipt_v0.1.json`
- `benchmarks/frozen/ttf_q_c2_result_receipt_v0.1.json`
- `benchmarks/frozen/ttf_q_b1_c1_information_ladder_v0.1.json`
- `benchmarks/frozen/ttf_q_bc_detectability_result_v0.1.json`
- `benchmarks/frozen/ttf_q_known_truth_decomposition_v0.1.json`
- `benchmarks/frozen/ttf_q_known_truth_detectability_v0.2.json` (after v0.2 result freeze)
- relational v0.3 remains `CLOSED_NO_EVALUABLE_TEST`.
