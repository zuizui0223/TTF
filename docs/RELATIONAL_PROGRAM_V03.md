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

For transparency, a five-test Holm calculation across S1-S3+B+C is reported only as a retrospective sensitivity label. Its first Holm threshold is 0.01 and it does not govern the primary v0.3 decision.

## Freshness

Study B retains the already acquired fresh-1000 external-data universe, but before pooled PCA or panel hashing it removes every species overlapping either prior confirmatory design universe. The exact response-blind audit found 25 overlaps with the S1 339-species design and 45 with the S2 500-species design, for **70 unique exclusions**. They are not backfilled.

Study C excludes the entire Study-B fresh-1000 in addition to all prior confirmatory design universes. Thus B and C cannot become species-overlapping because of a B gate failure or result.

## Study B

The environmental relation remains Schoener's D in a common two-axis whitened PCA of CHELSA V2.1 bio1, bio7, bio12 and bio15. Geography, class/order/family proximity, locality-count imbalance, and source/target fixed effects remain controls.

The qualification is now locked to the inferential alpha actually used by the family:

- p-value cutoff: **0.025**
- private-null Wilson 95% upper gate: **0.05**
- positive-control Wilson 95% lower power gate: **0.80**
- minimally meaningful standardized latent relation effect: **0.03**, unchanged

## Study C

C is frozen before B qualification. It uses CHELSA-TraCE21k centennial bioclimatic reconstructions from 21 ka BP to 0 BP. For bio01, bio07, bio12 and bio15, each occurrence receives a 21 ka minus 0 BP displacement vector. Species are represented by KDEs in a pooled, whitened two-axis PCA of those displacement vectors, and source-target historical similarity is Schoener's D.

The primary C model tests historical-displacement similarity beyond current CHELSA V2.1 environmental similarity, geographic opportunity, lineage controls, locality-count imbalance, and source/target fixed effects. Ice-sheet surface altitude is descriptive only and cannot become a rescue predictor.

## Stopping states

If B rejects, v0.3 closes and C is not opened. If B is null or NOT_EVALUABLE, C may proceed under its already frozen rule. After C, the program closes permanently.

Only two evaluable qualified nulls license the full two-relation null conclusion. A null plus a gate failure is reported as a closed partial-NOT_EVALUABLE program, not as two biological nulls.

Exact machine-readable control is in `docs/supporting/relational_program_v0.3.json`.
