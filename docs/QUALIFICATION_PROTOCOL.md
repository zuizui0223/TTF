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

The fixed-geometry Gate-I runner extends the same rule one level upward: because coordinates and the train/evaluation split are identical across semi-synthetic worlds, kNN graphs and the complete edge-integrated kernel projection are prepared once. Each world changes only synthetic trait values, turnover ranks, and species-bootstrap draws. A dedicated equivalence test requires this prepared path to reproduce the direct held-out-species pipeline.

## Gate G — adversarial type-I control

Mandatory null cells include:

- shared fraction = 0;
- amplitude > 0 at multiple strengths.

Every species may have a strong spatial transition, but transition locations are private.

The first prospective point-estimate criterion was:

- nominal alpha = 0.05;
- maximum rejection rate across zero-shared positive-amplitude cells <= 0.10.

The high-precision v0.2 criterion is stricter:

- 500 independently simulated worlds per mandatory cell;
- 1,999 held-out-species bootstrap resamples per world;
- the **two-sided Wilson 95% upper bound** for every zero-shared positive-amplitude rejection rate must be <= 0.10.

### Frozen v0.1 result — failed

Trait-permutation sharedness inference failed the original gate: maximum rejection reached **0.14** while fully shared moderate-signal power was 1.00. The failure is retained in `results/qualification_trait_permutation_pilot_v0.1.json`.

### Frozen v0.2 pilot — passed provisionally

Held-out-species bootstrap inference used 100 worlds per cell and 1,999 bootstrap resamples. Zero-shared rejection at amplitudes 0.5, 1, 2, 3 was **0.01, 0.07, 0.05, 0.07**; maximum = **0.07**. The pilot established that type I no longer increased with private-boundary amplitude, but 100 worlds were not sufficient for a narrow confidence bound.

### Frozen v0.2 high-precision result — passed

Workflow run `34107193204`, frozen head `cac260855e7846cb69e79ab2f79944132d721c11`, evaluated 500 worlds per key cell. Zero-shared rejection rates at amplitudes 0.5, 1, 2, 3 were:

- **0.040**, Wilson 95% = [0.0260, 0.0610];
- **0.050**, Wilson 95% = [0.0341, 0.0728];
- **0.040**, Wilson 95% = [0.0260, 0.0610];
- **0.054**, Wilson 95% = [0.0374, 0.0774].

The worst upper confidence bound was **0.07743**, below the prospectively frozen ceiling of 0.10. Gate G therefore **passes at high precision under the idealized calibration design**.

The complete result and workflow lineage are frozen in `results/qualification_heldout_species_precision_v0.2.json`.

## Gate H — positive-control power

Mandatory positive control:

- shared fraction = 1;
- prospectively defined moderate amplitude = 2.

The high-precision criterion requires the **two-sided Wilson 95% lower bound** on power to be >= 0.80.

In the same 500-world qualification, the fully shared amplitude-2 cell rejected in **489/500 = 0.978** worlds, Wilson 95% = **[0.9610, 0.9877]**. Gate H therefore **passes at high precision**.

A method that controls type I but has essentially no recovery of a moderate fully shared signal is uninformative, not validated.

## Gate I — empirical-geometry calibration

After idealized worlds pass, semi-synthetic calibration must be repeated on real sampling geometries. Trait values remain synthetic. This is the bridge between an abstractly functioning estimator and the identifiability available under empirical sampling.

Gate I is deliberately split into two stages because they answer different failure modes.

### Gate I-A — independent non-flower geometry stress

Purpose: show that qualification is not an artefact of the flower-colour sampling geometry used by the motivating research programme.

The first frozen source manifest uses 16 widespread European bird and mammal taxa. The acquisition path:

- uses coordinate-bearing GBIF occurrences in a prospectively declared European bounding box;
- rejects unresolved taxon matches and records with geospatial issues;
- removes only literal duplicate coordinate pairs;
- applies a deterministic SHA-256 source-key priority when a species exceeds the frozen cap;
- converts latitude/longitude to 3-D Earth-centred kilometre coordinates so Euclidean graph distances are meaningful over the continental domain;
- freezes retained GBIF occurrence-key fingerprints and the resulting TTF geometry fingerprint;
- selects the smoothing bandwidth from geometry only, before synthetic trait outcomes are generated.

The small external workflow pilot validates acquisition and execution only. It does **not** count as a high-precision qualification.

### Gate I-B — intended empirical geometry

Purpose: establish identifiability on the exact sampling frame intended for the eventual TTF empirical claim.

The following must be retained exactly from that intended frame:

- species counts and identities;
- coordinates;
- graph geometry;
- record counts after prospectively declared QC/capping;
- train/evaluation constraints;
- opportunity structure;
- dependence blocks when required.

A successful I-A result cannot substitute for I-B. Conversely, I-B should not be used as the only cross-domain stress test when an independent geometry is available.

### Gate-I implementation contract

The fixed-geometry engine is implemented in `ttf.geometry` and `run_geometry_calibration`.

Every Gate-I run must satisfy all of the following:

- empirical trait values are not passed into the simulator; the input object stores only species identity, coordinates and optional dependence blocks;
- raw coordinates, species identities and record counts are copied unchanged into every synthetic world;
- the train/evaluation split is frozen prospectively and reused across worlds unless the intended deployment protocol explicitly specifies another split schedule;
- synthetic shared transitions use one common geographic hyperplane, while zero-shared controls use species-private transition hyperplanes;
- affine standardization is used only to define synthetic transition strength in unit-free coordinates; graph construction, distances, bandwidth and boundary-field estimation continue to use the original empirical coordinates;
- the exact sampling frame is frozen with a SHA-256 geometry fingerprint and per-species record-count ledger;
- kNN graphs and the edge-integrated kernel projection are geometry-only objects prepared once and reused across worlds;
- the prepared fast path must reproduce direct inference numerically;
- the primary inference remains the held-out-species bootstrap; trait permutation is not substituted for Gate-I sharedness inference.

The claim-bearing I-A and I-B qualification runs use the same high-precision confidence-bound criteria as Gates G/H: 500 worlds per mandatory cell and 1,999 bootstrap resamples, with every zero-shared positive-amplitude Wilson 95% upper bound <= 0.10 and the fully shared moderate-amplitude Wilson 95% lower bound >= 0.80.

**Implementation status:** fixed-geometry engine, fast path, CI smoke tests, external animal manifest, and external acquisition/pilot workflow exist. **Qualification status:** neither I-A nor I-B has yet passed a 500-world high-precision run. Gate I therefore remains open.

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

**Passed at high precision on idealized circular worlds:**

- species isolation and estimator invariants;
- species-disjoint edge-level transfer;
- held-out-species inferential-unit tests;
- strong-private-boundary type-I gate, including Wilson uncertainty;
- fully shared moderate-signal power gate, including Wilson uncertainty.

**Implemented but not yet empirically qualified:**

- fixed empirical-geometry semi-synthetic engine and lineage fingerprinting;
- Gate-I geometry-only fast path;
- Gate-I-A external non-flower acquisition/pilot path.

Still downstream and explicitly unqualified:

- a completed high-precision Gate-I-A external geometry qualification;
- a completed high-precision Gate-I-B intended-geometry qualification;
- genetic-distance isolation-by-distance residualization;
- environmental/resistance feature builders;
- full observation-bias stress suite;
- a cross-domain empirical detection floor.

The method-core sharedness inference is therefore qualified for its **idealized synthetic scope**, not yet for an arbitrary empirical dataset.
