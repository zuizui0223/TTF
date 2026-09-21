# Relational TTF v0.3 — finite prospective program

Status: **frozen before Study B formal synthetic qualification**.

## What changed

The ecological-predictor search now has a finite end. Three completed predictor tests remain immutable prior evidence: butterfly trait similarity (p=0.89), host-resource geography (p=0.698), and relational host-resource similarity (p=0.6626). S2 and S3 are not described as independent species confirmations; after strict parser repair, the 420 S3 survivors are contained within the 431 S2 survivors.

The new prospective v0.3 family has exactly two confirmatory slots:

- **B — environmental niche similarity**, one-sided alpha 0.025.
- **C — historical climate displacement similarity**, one-sided alpha 0.025.

No third ecological predictor may be added. Gate failure is NOT_EVALUABLE, consumes the slot, and does not transfer its alpha to the other slot.

## Multiplicity

The primary v0.3 family uses fixed Bonferroni allocation: 0.025 + 0.025 = 0.05. This is a new prospective family definition; it does **not** claim that prior non-significant tests had no Type-I error opportunity.

For transparency, a five-test Holm calculation across S1-S3+B+C is reported only as a retrospective sensitivity label **when both B and C empirical p-values exist**. Its first Holm threshold is 0.01 and it does not govern the primary v0.3 decision. If B rejects and the stopping rule leaves C unopened, no five-test Holm result is reported and no p-value is imputed for C.

## Freshness

Study B retains the already acquired fresh-1000 external-data universe, but before pooled PCA or panel hashing it removes every species overlapping either prior confirmatory design universe. The exact response-blind audit found 25 overlaps with the S1 339-species design and 45 with the S2 500-species design, for **70 unique exclusions**. They are not backfilled.

Study C excludes the entire Study-B fresh-1000 in addition to all prior confirmatory design universes. The exclusions are not hand-curated: S1 is regenerated and required to match its frozen design SHA, S2 is regenerated from the frozen archive plus WCVP/HOSTS sidecars and required to match its selected-500 digest, and the B candidate CSV must match its frozen SHA. Their verified set union is excluded before the C hash ranking. Thus B and C cannot become species-overlapping because of a B gate failure or result.

## Study B

The environmental relation remains Schoener's D in a common two-axis whitened PCA of CHELSA V2.1 bio1, bio7, bio12 and bio15. Geography, class/order/family proximity, locality-count imbalance, and source/target fixed effects remain controls.

The machine-readable B chain is `relational_environment_relation_rule_v0.3.json` → `relational_environment_opportunity_rule_v0.2.json` → `relational_environment_qualification_rule_v0.2.json`. The opportunity rule changes no geography threshold; v0.2 exists only to bind the unchanged 500-km / 0.50-coverage / five-source contract to the exact post-freshness v0.3 panel.

The qualification is now locked to the inferential alpha actually used by the family:

- p-value cutoff: **0.025**
- private-null Wilson 95% upper gate: **0.05**
- positive-control Wilson 95% lower power gate: **0.80**
- minimally meaningful standardized latent relation effect: **0.03**, unchanged

## Study C

C is frozen before B qualification. It uses CHELSA-TraCE21k centennial bioclimatic reconstructions from 21 ka BP to 0 BP. For bio01, bio07, bio12 and bio15, each occurrence receives a 21 ka minus 0 BP displacement vector. Species are represented by KDEs in a pooled, whitened two-axis PCA of those displacement vectors, and source-target historical similarity is Schoener's D. The frozen logical assets are the eight V1.0 variable/time combinations at indices -190 and 20; the exact provider URLs and SHA256 values are recorded at the first response-blind asset resolution, with no scientific asset substitution allowed.

The C candidate scan inherits the exact response-blind geometry thresholds: Animalia excluding Aves/Chiroptera, COI-family marker, at least 12 unique localities, neighbor fraction 0.15, and at least five endpoint-disjoint IBD training edges. It retains the first 1000 eligible species in the frozen history hash order, or all 500-999 if fewer than 1000 remain; fewer than 500 closes C as NOT_EVALUABLE.

The primary C model tests historical-displacement similarity beyond current CHELSA V2.1 environmental similarity, geographic opportunity, lineage controls, locality-count imbalance, and source/target fixed effects. Its separate `relational_historical_climate_qualification_rule_v0.1.json` freezes alpha 0.025, 1000 worlds per cell, the A=0.5/1/2/3 private-null family, source/target intercept SD 0.18, dyad-noise SD 0.32, beta_hist floor 0.03, Wilson Type-I upper gate 0.05, and power lower gate 0.80. Ice-sheet surface altitude is descriptive only and cannot become a rescue predictor.

## Stopping states

If B rejects, v0.3 closes and C is not opened. If B is null or NOT_EVALUABLE, C may proceed under its already frozen rule. After C, the program closes permanently.

Only two evaluable qualified nulls license the full two-relation null conclusion. A null plus a gate failure is reported as a closed partial-NOT_EVALUABLE program, not as two biological nulls.

Exact machine-readable control is in `docs/supporting/relational_program_v0.3.json`.
