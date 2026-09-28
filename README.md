# TTF — Transferability Qualification and Abstention

TTF is a methods repository for a single question:

> When may structure learned in some biological units be claimed to transfer to entirely unseen units, and when should the analysis abstain instead?

The repository develops prospective qualification, held-out evaluation, finite candidate families, frozen stopping rules, and explicit `NOT_EVALUABLE` endpoints for cross-unit transfer claims. Biological systems are used as test beds for the methodology; they are not separate ecological projects within this repository.

## Core principle

TTF separates four evidence layers that are often conflated:

```text
method validity on declared synthetic worlds
        !=
deployment-geometry detectability
        !=
observational admissibility
        !=
empirical biological conclusion
```

A method may be statistically valid yet unusable for a particular dataset. A dataset may be observationally admissible yet show no cross-unit transfer. A failed qualification is therefore not an ecological null result: it is a stopping decision.

## TTF-Q: evaluable envelopes for dyadic ecological predictors

TTF-Q is the active response-blind method-development route for pairwise ecological predictors. It does **not** reopen the closed relational v0.3 B/C confirmatory family. Instead, it asks what can be established about a source-target relation before any biological response is opened.

For each focal relation, TTF-Q decomposes:

- variation removed by exact source and target identity;
- additional redundancy with declared ecological and sampling controls;
- information attributable to geographic co-opportunity;
- concentration of residual relational signal across dyads, sources, and targets;
- the effect-size × nuisance-heterogeneity region in which the intended dyadic estimator is detectably informative.

The output is an **evaluable envelope**, not a binary PASS/FAIL label.

The current response-blind climate case study shows two distinct information layers on a shared dyad geometry. Present-climate similarity retains non-redundant pairwise information after endpoint identity, lineage, sampling imbalance, and geography. Historical climate-displacement similarity adds a second layer beyond present climate; after adding continuous geographic co-opportunity, about **82–83%** of that previously unique historical information remains across the two frozen panels.

These are ecological-predictor information results, not tests of a genetic-transferability effect.

Primary TTF-Q materials:

- `docs/TTF_Q_V01.md` — method definition;
- `docs/TTF_Q_MANUSCRIPT_SPINE_V01.md` — manuscript logic and claim boundaries;
- `benchmarks/frozen/ttf_q_c1_result_receipt_v0.1.json` — history beyond current climate;
- `benchmarks/frozen/ttf_q_c2_result_receipt_v0.1.json` — geography decomposition;
- `benchmarks/frozen/ttf_q_b1_c1_information_ladder_v0.1.json` — present-to-historical information ladder;
- `benchmarks/frozen/ttf_q_bc_detectability_result_v0.1.json` — detectability envelopes.

## Primary estimand

For a held-out unit (s), TTF learns a predictor from training units only and evaluates it on data from units excluded from fitting.

For the original turnover-field formulation,

```text
C_s = Spearman(predicted transition structure, observed transition structure)
T   = mean_s C_s over held-out units
```

The exact response and nuisance adjustment vary by application, but the cross-unit logic does not:

1. define the estimand before opening held-out outcomes;
2. freeze admissibility and qualification rules;
3. keep training and evaluation units disjoint;
4. calibrate error and power on declared worlds;
5. stop when the design is not evaluable;
6. do not rescue a failed frozen design by outcome-dependent retuning.

## Current flagship empirical route

The completed genetic route asks whether within-species genetic differentiation learned from training species predicts differentiation in unseen species beyond the declared isolation-by-distance adjustment.

The terminal fresh panel contains 211 eligible species after response-blind filtering, with an inherited 103/108 training/evaluation split. Synthetic Gate-D and self-detectability qualification passed before the empirical response was opened.

The frozen empirical endpoint is:

- cross-species post-IBD transfer statistic: (T = 0.03257), (p = 0.7393);
- calibrated within-species self-detectability test: (p = 0.00899);
- terminal decision: `LINEAGE_CONDITIONED_SPATIAL_STRUCTURE_WITHIN_TESTED_DOMAIN`.

This is a qualified non-detection of cross-species transfer under the declared estimand, alongside detectable within-species structure relative to the frozen structural null. The raw self statistic is negative, so it must not be described as a positive raw correlation.

Primary frozen handoff:
`benchmarks/frozen/genetic_phylogatr_phase4_empirical_handoff_v0.1.json`

Current manuscript:
`manuscript/genetic_ttf_flagship_v0.2.md`

## Why abstention is part of the method

TTF treats `NOT_EVALUABLE` as a valid terminal result rather than an invitation to tune until a signal appears.

Examples already frozen in the repository include:

- a development design that failed its private-structure error-control bound;
- an observation-support gate in which too few species retained the required number of evaluable records;
- heterogeneous-geometry relation-transfer candidates for which no prospectively declared transport passed all development margins;
- relational candidates that were closed when required external data could not support the predeclared estimand.

These records are retained because the ability to stop is part of the methodological claim.

## Finite-family discipline

Later candidate families are explicitly finite. A candidate is admitted, qualified, or rejected under rules frozen before the corresponding response is opened. Once a family is exhausted, the program closes rather than generating post-outcome replacements under the same confirmatory label.

This repository therefore preserves failed gates and closure receipts alongside successful qualification. They document the selection process required to interpret the final empirical test.

## Repository boundary

Ecological analyses whose scientific contribution is the biology rather than the transfer methodology are being moved out of TTF.

The butterfly specialization / host-resource / climate analysis is migrating to:

`https://github.com/zuizui0223/chocho`

The migration is pinned to TTF commit:

`1a112334cca2f2ef5e234c3ae1fc1a80b8266956`

The butterfly files remain on TTF `main` until the independent repository passes its paper-specific test suite. They will then be removed from the methods-only branch without rewriting TTF history.

## Repository map

- `src/ttf/` — estimators, nulls, qualification gates, simulation and frozen execution logic.
- `tests/` — implementation and frozen-result tests.
- `docs/supporting/` — prospectively frozen rules, amendments, and closure documents.
- `benchmarks/frozen/` — immutable qualification and empirical receipts.
- `manuscript/` — methods manuscript material and generated handoffs.
- `scripts/` — reproducible builders, qualification runs, audits, and renderers.

Exploratory biological analyses are not part of the long-term methods-only surface and are being separated into paper-specific repositories.

## Reproducibility contract

A result intended to support a manuscript claim should be traceable through:

```text
frozen rule
  -> source/data identity
  -> executable implementation
  -> test or qualification gate
  -> immutable result receipt
  -> manuscript claim
```

The repository should make it possible to determine not only how a reported result was obtained, but also which alternatives were eligible before that result was known and why the analysis stopped where it did.

## Development status

The methods-only refactor is in progress. No historical frozen result is being rewritten. The active cleanup removes ecological paper surfaces only after their destination repository reproduces them independently.
