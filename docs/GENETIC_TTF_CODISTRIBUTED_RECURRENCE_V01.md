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

## Fixed-dyad geometry-null implementation status

The next inferential layer is now frozen in `docs/supporting/genetic_codistributed_recurrence_geometry_null_rule_v0.3.json` and implemented in `src/ttf/codistributed_geometry_null.py`.

For each directed dyad, the geometry-null center is

`mu0_st = mean_A mean_r T_st(private world A, replicate r)`,

with equal weight assigned to each predeclared private amplitude `A in {0.5, 1, 2, 3}`. The primary response is then `E_st = T_st - mu0_st`.

This subtraction is not used as a parametric null assumption. The final one-sided decision remains calibrated against independent private-world beta_G reference distributions at each amplitude, with the largest component Monte Carlo p-value taken as the primary p-value. Thus an incorrect guess about the empirical private amplitude cannot be selected after response opening.

Screen, formal-reference, formal-evaluation and shared-positive worlds use disjoint seed namespaces fixed before any empirical T_st is opened. The actual source-to-target operator is the inherited endpoint-safe post-IBD, 500-km single-source Gaussian TTF operator; a batch implementation is tested against the prior scalar implementation.

Execution is currently blocked only by absence of the exact 274,988,692-byte phylogatR archive (SHA-256 `5a0fd9ac25893c749d14186fbcce4a46b99163c9d810b36e40eebce7bece61a5`). The retained compact midpoint geometry is sufficient for the response-blind G_st / TTF-Q audits but not for reconstructing endpoint-safe graph nodes required by the actual fixed-dyad T_st qualification. No genetic response has been opened.

## Comparative-phylogeographic novelty

The study is intentionally anchored in the classical **codistributed-species** problem rather than in an animal-wide claim about universal barriers.

Previous comparative phylogeography has mainly asked whether already-observed taxa show spatially concordant breaks, divergence histories, or barrier associations. Examples include congruence tests among nine Californian codistributed taxa (Lapointe & Rissler 2005), explicit tests of microhabitat-mediated discordance between two codistributed montane sedges (Massatti & Knowles 2014), and the recent synthesis of mapped phylogeographic breaks across 229 North American mammal species (Jensen et al. 2024). Papadopoulou & Knowles (2016) argued that such discordance should itself become predictable from organismal biology.

The present prospective study asks a different question:

> **Does the amount of sampled geography actually shared by two lineages predict how much post-IBD spatial genetic information transfers from one lineage to the other?**

This distinction matters. `G_st` is measured before genetic response, continuously for every fixed source-target dyad, and the target lineage remains absent from source-field fitting. The response is predictive reuse of place information, not retrospective overlap of manually identified breaks and not proximity to a pre-labelled mountain, river, or ecotone.

The closest large-scale contrast is Jensen et al. (2024): their analysis asks where literature-derived breaks in 229 North American mammals fall relative to candidate geographic barriers and whether species traits explain variation in barrier effects. Here no break catalogue or named barrier layer enters the primary estimand. Instead, geographic co-opportunity itself is asked to predict out-of-lineage transfer skill across a fresh prospective insect panel.

This yields a sharper biological distinction:

- **place component:** greater geographic co-opportunity raises excess source-to-target recurrence;
- **lineage-conditioned component:** even when two lineages share sampled geography, their post-IBD spatial differentiation remains non-reusable.

A positive `beta_G` would therefore support a reusable place component without claiming that a specific physical barrier caused it. A qualified null would show that simple lack of shared geographic opportunity is insufficient to explain the broad 211-species non-transfer result.

### References added for positioning

- Lapointe, F.-J. & Rissler, L. J. (2005). Congruence, consensus, and the comparative phylogeography of codistributed species in California. *The American Naturalist* 166:290–299. doi:10.1086/431283.
- Massatti, R. & Knowles, L. L. (2014). Microhabitat differences impact phylogeographic concordance of codistributed species: genomic evidence in montane sedges from the Rocky Mountains. *Evolution* 68:2833–2846. doi:10.1111/evo.12491.
- Papadopoulou, A. & Knowles, L. L. (2016). Toward a paradigm shift in comparative phylogeography driven by trait-based hypotheses. *PNAS* 113:8018–8024. doi:10.1073/pnas.1601069113.
- Jensen, A. J. et al. (2024). Geographic barriers but not life history traits shape the phylogeography of North American mammals. *Global Ecology and Biogeography* 33:e13875. doi:10.1111/geb.13875.

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
