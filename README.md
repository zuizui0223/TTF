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

The ecological butterfly paper has been separated from TTF into:

`https://github.com/zuizui0223/chocho`

The split is pinned to TTF commit:

`1a112334cca2f2ef5e234c3ae1fc1a80b8266956`

The destination repository independently passes its paper-specific offline test suite, and the butterfly/Lepidoptera application surface has been removed from current TTF `main`. TTF history is unchanged, so the pre-split files remain inspectable through the frozen source commit.

## Repository map

- `src/ttf/` — estimators, nulls, qualification gates, simulation and frozen execution logic.
- `tests/` — implementation and frozen-result tests.
- `docs/supporting/` — prospectively frozen rules, amendments, and closure documents.
- `benchmarks/frozen/` — immutable qualification and empirical receipts.
- `manuscript/` — methods manuscript material and generated handoffs.
- `scripts/` — reproducible builders, qualification runs, audits, and renderers.
- `.github/workflows/` — current CI and deliberately supported reproducibility entrypoints.
- `provenance/workflows/` — completed, superseded, diagnostic, and recovery workflow definitions preserved byte-for-byte as historical provenance.

The workflow classification is machine-readable in `provenance/workflow_surface_classification_v0.1.json`. Exploratory biological papers live in paper-specific repositories rather than the active TTF surface.

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

The methods-only repository split is complete. No historical frozen result was rewritten. The current cleanup phase reduces active repository surface area while preserving historical execution definitions and frozen receipts under explicit provenance paths.

Repository-wide CI is `.github/workflows/tests.yml`. Historical workflow moves preserve the original Git blob identity and therefore do not alter the recorded execution definitions.
