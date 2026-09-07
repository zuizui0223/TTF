# TTF Qualification Protocol v0.2

A TTF empirical result is not claim-ready merely because code executes or a p-value is small. Qualification is prospective and the inference target must be stated before empirical outcomes are inspected.

## Gate A — species isolation

**Requirement:** graph construction is species-local by construction. No graph utility may silently connect records from different species.

**Automated test:** species-local graph invariant.

## Gate B — rank and scheduling invariants

**Requirements:**

- within-species ranks are invariant under strictly monotone rescaling of dissimilarity;
- no species is duplicated inside one realization;
- long-run inclusion counts differ by at most one when balanced recurrence schedules are used.

## Gate C — species-disjoint prediction

**Requirements:**

- training and evaluation species are disjoint;
- the training field never uses evaluation trait values;
- each evaluation species is scored from edge-level predictions;
- the primary statistic equals the macro-average of finite held-out species scores.

Sharedness is defined by this transfer step, not by a pooled hotspot or map-map similarity.

## Gate D — primary inferential-unit integrity

The v0.2 sharedness null is conditional on the learned training field:

\[
H_0:E[C_s\mid B_{\rm train}]\le0.
\]

**Requirements:**

- held-out species, not records or edges, are the resampling units;
- coordinates, trait values, graph geometry, and private within-species spatial structure are untouched during the primary sharedness bootstrap;
- centered bootstrap resampling imposes a zero mean without erasing between-species heterogeneity;
- the bootstrap implementation is deterministic under a fixed seed.

This is a superpopulation inference about unseen species conditional on the trained field.

## Gate E — diagnostic permutation invariance

Within-species trait permutation remains available for the distinct question of trait-location exchangeability.

When that diagnostic is used:

- coordinates and graph geometry must remain unchanged;
- only allowable trait labels are shuffled;
- dependence blocks must be respected when present;
- pairwise distance-matrix cells must never be shuffled independently.

A diagnostic trait-permutation p-value must not be reported as the primary sharedness p-value.

## Gate F — exact fast-path equivalence

Geometry-only kernel precomputation is accepted only if it reproduces direct edge-integrated boundary prediction to numerical precision.

This prevents estimator drift from applying opportunity correction after averaging along an edge rather than at each integration point.

## Gate G — adversarial type-I control

Mandatory null cells include:

- shared fraction = 0;
- amplitude > 0 at multiple strengths.

Every species may have a strong spatial transition, but transition locations are private.

Prospective provisional criterion:

- nominal alpha = 0.05;
- maximum rejection rate across zero-shared positive-amplitude cells <= 0.10.

### Frozen v0.1 result

Trait-permutation sharedness inference failed this gate: maximum rejection reached **0.14** while fully shared moderate-signal power was 1.00. The failure is retained in `results/qualification_trait_permutation_pilot_v0.1.json`.

### Frozen v0.2 pilot result

Held-out-species bootstrap inference used 100 worlds per cell and 1,999 bootstrap resamples. Zero-shared rejection at amplitudes 0.5, 1, 2, 3 was **0.01, 0.07, 0.05, 0.07**; maximum = **0.07**, passing the provisional point-estimate ceiling. The result is retained in `results/qualification_heldout_species_bootstrap_v0.2.json`.

This is still a pilot: Wilson 95% upper bounds for a 0.07 rate with 100 worlds reach about 0.137. Higher-replicate type-I calibration is therefore required before calling the method publication-grade.

## Gate H — positive-control power

Mandatory positive control:

- shared fraction = 1;
- prospectively defined moderate amplitude.

Provisional target: power >= 0.80.

The v0.2 pilot gives **0.96 power at amplitude 2**, passing this gate.

A method that controls type I but has essentially no recovery of a moderate fully shared signal is uninformative, not validated.

## Gate I — actual-geometry calibration

After idealized worlds pass, repeat semi-synthetic calibration using the intended empirical:

- species counts;
- coordinates;
- graph geometry;
- record counts;
- train/evaluation constraints;
- opportunity structure.

Trait values remain synthetic. This is the bridge between an abstractly functioning estimator and the identifiability available in a real opportunistic dataset.

## Gate J — predictor attribution

Geographic, environmental, and resistance predictor spaces must be compared on identical held-out species. Report paired held-out increments such as `T_GE - T_G` rather than winner-take-all labels.

For traits with intrinsic isolation-by-distance, nuisance IBD must be cross-fit before shared residual structure is evaluated.

## Gate K — observation-bias stress tests

Before a broad geographic claim:

- freeze availability surfaces;
- cap dominant observers or sources;
- delete high-effort cells;
- run location-blind measurement where feasible;
- check whether negative-control availability predicts the learned boundary field.

## Current v0.2 status

Passed provisionally on idealized circular worlds:

- species isolation and estimator invariants;
- species-disjoint edge-level transfer;
- held-out-species inferential-unit tests;
- strong-private-boundary point-estimate type-I gate;
- fully shared moderate-signal power gate.

Not yet qualified:

- publication-grade type-I precision;
- actual opportunistic sampling geometry;
- genetic-distance IBD residualization;
- environmental/resistance feature builders;
- full observation-bias stress suite;
- a publication-grade detection floor.

These remain explicit downstream gates rather than being hidden inside an empirical analysis.
