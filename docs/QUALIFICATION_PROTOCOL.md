# TTF Qualification Protocol v0.1

A TTF empirical result is not claim-ready merely because code executes or a permutation p-value is small. Qualification is prospective.

## Gate A — species isolation

**Requirement:** graph construction is species-local by construction. No graph utility accepts pooled species labels and silently connects them.

**Automated test:** `test_graphs_are_species_local_by_construction`.

## Gate B — rank and scheduling invariants

**Requirements:**

- within-species ranks are invariant under strictly monotone rescaling of dissimilarity;
- no species is duplicated inside one realization;
- long-run inclusion counts differ by at most one.

**Automated tests:** rank-scale and balanced-schedule tests.

## Gate C — null geometry invariance

**Requirement:** trait permutation may change turnover values but not coordinates, graph nodes, edge positions, environmental geometry, or the species split.

**Automated test:** coordinates and edge geometry are byte-for-byte unchanged after permutation.

When dependence blocks are supplied, shuffling is restricted to those blocks.

## Gate D — exact fast-path equivalence

Large permutation counts require geometry-only kernel precomputation. Optimization is accepted only if it reproduces direct edge-integrated boundary prediction.

**Automated test:** prepared geometry and direct scoring agree to numerical precision.

This gate prevents a subtle estimator drift where opportunity correction is applied after averaging along an edge instead of at each integration point.

## Gate E — adversarial type-I control

Mandatory cells include:

- shared fraction = 0;
- amplitude > 0 at multiple strengths.

Interpretation: every species may have a strong spatial transition, but transition locations are private.

Prospective provisional criterion:

- nominal alpha = 0.05;
- maximum rejection rate across zero-shared positive-amplitude cells <= 0.10.

Before publication, this should be replaced or supplemented by a binomial confidence bound sized to the number of calibration replicates.

## Gate F — positive-control power

Mandatory positive control:

- shared fraction = 1;
- prospectively defined moderate amplitude.

Provisional target:

- power >= 0.80.

A method that controls type I but has essentially no recovery of a moderate fully shared signal is **uninformative**, not validated.

## Gate G — actual-geometry calibration

After idealized circular worlds pass, repeat semi-synthetic calibration using the empirical:

- species counts;
- coordinates;
- graph geometry;
- record counts;
- balanced scheduling constraints.

Trait values remain synthetic. This is the bridge between an abstractly valid estimator and the identifiability available in a real opportunistic dataset.

## Gate H — predictor attribution

Geographic, environmental, and resistance predictor spaces must be compared on identical held-out species and identical null replicates. Report paired increments such as `T_GE - T_G`.

## Gate I — observation-bias stress tests

Before a broad geographic claim:

- freeze availability surfaces;
- cap dominant observers or sources;
- delete high-effort cells;
- run location-blind measurement where feasible;
- check whether negative-control availability predicts the learned boundary field.

## Current v0.1 status

The repository contains the core estimator, full permutation null, idealized calibration engine, predictor-space transfer interface, and invariant tests.

Not yet qualified in v0.1:

- actual RGFCA-like sampling geometry;
- genetic-distance IBD residualization;
- environmental/resistance feature builders;
- full observation-bias stress suite;
- a publication-grade detection floor.

Those are deliberately downstream of the type-I/power gate rather than prerequisites hidden inside an empirical analysis.
