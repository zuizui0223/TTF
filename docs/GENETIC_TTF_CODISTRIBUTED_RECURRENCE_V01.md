# Codistributed recurrence in fresh Insecta v0.1

Status: **prospective response-blind design frozen before genetic response opening**

## Question

> Does sharing more of the same landscape make post-IBD phylogeographic structure more recurrent across lineages?

This directly tests the strongest alternative explanation for the completed 211-species result: perhaps cross-lineage recurrence was weak simply because many lineages did not experience enough of the same places.

## Fresh domain

Start from the frozen response-blind fresh-1000 geometry set from the exact authenticated phylogatR archive.

- original fresh 250 species already excluded;
- prior conditional 500 species already excluded;
- additionally exclude the 70 fresh-1000 species overlapping frozen S1/S2 response-opening universes;
- retain Insecta only.

This leaves **762 response-blind insect species**. No retained species has contributed genetic response to the completed TTF analyses used to motivate this hypothesis.

## Deterministic panels

Rank species by:

`SHA256("place-recurrence-insecta-fresh-v0.1|<archive-sha256>|<species>")`

First 381 = development; remaining 381 = confirmatory.

Separate role hashes assign 190 source and 191 target species within each panel.

Development nucleotide identity remains closed permanently. Confirmatory identity may be opened once only after all gates pass.

## Primary relation

For fixed source s and target t, use graph-edge midpoint support at the inherited 500-km scale.

- target-to-source coverage = fraction of target midpoints within 500 km of a source midpoint;
- source-to-target coverage = reverse quantity;
- primary continuous relation:

`G_st = min(target_to_source, source_to_target)`.

All 36,290 directed dyads per panel are retained. No high/low threshold is used in the primary model.

## Response-blind feasibility

For interpretation only, using `G_st >= 0.50` and requiring at least five sources per target:

- development: **162/191 targets**, median 32 supported sources;
- confirmatory: **158/191 targets**, median 25 supported sources.

The supported target set spans multiple insect orders rather than one taxonomic subgroup.

## TTF-Q information audit

Nested controls:

1. exact source and target identity;
2. absolute log edge-count ratio;
3. absolute log locality-count ratio;
4. same-order indicator;
5. log1p source-target centroid distance / 500 km.

### Development

- endpoint survival: **0.8636**;
- total survival after sampling + lineage controls: **0.8484**;
- final total unique fraction after centroid distance: **0.2149**;
- effective sources: **168.1/190**;
- effective targets: **156.8/191**;
- max source signal share: **1.19%**;
- max target signal share: **1.97%**.

### Confirmatory

- endpoint survival: **0.8585**;
- total survival after sampling + lineage controls: **0.8449**;
- final total unique fraction after centroid distance: **0.2054**;
- effective sources: **162.9/190**;
- effective targets: **158.9/191**;
- max source signal share: **1.47%**;
- max target signal share: **1.76%**.

Thus geographic co-opportunity is not merely endpoint identity, sampling imbalance, same-order membership, or centroid proximity. About 20–21% of its total variation remains uniquely relational after all declared controls in both panels, with broad endpoint support.

This is a predictor-information result only. No genetic association has been opened.

## Future response and estimand

For every fixed source-target dyad:

- `T_st` = post-IBD source-to-target transfer skill under one fixed operator;
- `mu0_st` = expected transfer under the frozen geometry-specific private-null family;
- `E_st = T_st - mu0_st`.

Primary model:

`E_st = source FE + target FE + beta_G G_st + declared controls + error_st`.

Directional hypothesis: **beta_G > 0**.

A positive result would mean that lineages sharing more of the same landscape show greater excess recurrence of spatial genetic structure than expected from lineage-private spatial worlds.

A calibrated null, together with qualified within-species structure, would show that lack of shared geographic opportunity is not sufficient to explain the broad 211-species non-transfer result.

## Required opening gates

Before confirmatory nucleotide identity:

1. implement geometry-null centering for the fixed dyad operator;
2. qualify two-way source/target inference on development synthetic worlds;
3. require Type-I control over private amplitudes A = 0.5, 1, 2, 3;
4. require power for a frozen shared-place positive control;
5. map a calibration-qualified detectable beta_G envelope;
6. open character-validity masks only for confirmatory species;
7. rerun the exact TTF-Q information and synthetic gates on survivor geometry.

Any failure yields `NOT_EVALUABLE_CODISTRIBUTED_RECURRENCE`.

No radius, threshold, subgroup, alternate overlap metric, or source-pool rule may be chosen after empirical response opening.

## Claim boundary

The response is recurrent place-specific differentiation, not a direct map of named physical barriers.

A positive result would support reusable geography among codistributed lineages, but would not identify whether the mechanism is a mountain, river, climatic transition, refugial boundary, secondary contact, or another process.
