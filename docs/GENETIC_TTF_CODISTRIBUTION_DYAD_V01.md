# Genetic TTF codistribution dyad successor v0.1

Status: **new response-blind successor design; no new empirical genetic response authorized**

Parent biological result: completed 211-species unconditional genetic TTF.

Parent interpretation: spatial genetic structure is detectable within species, but its post-IBD geographic configuration was not detectably reusable across unseen lineages under the tested broad multispecies field.

## 1. Scientific question

The unconditional result averages over lineages with very different geographic support. The next question is therefore narrower and more directly comparative-phylogeographic:

> **Among source and target lineages that share more of the same landscape, is post-IBD spatial genetic structure more recurrent than expected from lineage-private spatial structure alone?**

This asks whether **codistribution creates geographic recurrence**.

It does not ask whether a preselected taxonomic group has a positive transfer score, and it does not reopen any species-level outcome from the completed 211-species result.

## 2. Why the source-pool approach is closed

An earlier conditional program changed the training source pool separately for each target.

- v0.2 compared geographically supported same-order and different-order source pools.
- The development geometry qualified, but the exact confirmatory character-mask survivor geometry failed the private-null Type-I gate before nucleotide identity was opened.
- The response-blind diagnosis showed that same-order and different-order pools differed materially in geographic overlap, centroid distance, edge count and locality count.
- v0.3 prospectively replaced the zero-centered bootstrap with a least-favourable private-null envelope and matched same-order and different-order sources on response-blind geometry.

The v0.3 development screen is now numerically reproduced on the exact frozen development geometry under its already-frozen seed namespaces:

- private A0.5: 0/50 rejections;
- private A1: 1/50;
- private A2: 1/50;
- private A3: 0/50;
- same-order A2: 37/50;
- frozen development-screen criteria: private cells <=5/50 and positive cell >=35/50;
- development result: PASS.

A semantic reconstruction of the 214-species confirmatory survivor geometry then produced:

- private A0.5: 0/50;
- private A1: 0/50;
- private A2: 1/50;
- private A3: 0/50;
- same-order A2: 15/50.

Thus geometry matching appears to repair the v0.2 Type-I pathology but does not recover the predeclared positive-control detectability target on the reconstructed survivor geometry.

**Audit boundary:** the available survivor CSV is a semantic reconstruction, not the byte-identical canonical survivor artifact bound by the earlier receipt. Therefore 15/50 is a development diagnostic, not a formal frozen qualification result. It is sufficient to motivate not spending a fresh empirical opening on this source-pool estimator, but it is not reported as an empirical or canonical qualification claim.

Under the frozen v0.3 formalization rule, thresholds, matching features and positive-control strength must not be tuned against these worlds. The same-order source-pool route is therefore not advanced.

## 3. New principle: keep the response operator fixed

Changing the source pool confounds the biological relation with the geometry of the prediction operator.

The successor instead retains a fixed source-to-target transfer operator for every eligible directed dyad.

For source species s and target species t:

- T_st = post-IBD transfer skill from source s to target t under one frozen operator;
- G_st = response-blind geographic co-opportunity between source s and target t;
- P_st = predeclared lineage proximity control, if used;
- X_st = other geometry/sampling controls declared before response opening.

All eligible dyads remain in the response table. Geography is a predictor, not a rule for choosing which source response is allowed to enter the field.

## 4. Primary geographic relation

The preferred primary G_st is a continuous, symmetric edge-midpoint co-opportunity measure at the inherited 500-km scale.

For each directed pair, compute:

- forward coverage: fraction of target graph-edge midpoints with a source midpoint within 500 km;
- reverse coverage: fraction of source graph-edge midpoints with a target midpoint within 500 km.

Primary symmetric co-opportunity:

`G_st = min(forward_coverage, reverse_coverage)`.

Rationale:

- it requires both lineages to occupy spatially comparable support;
- it does not reward a tiny source range nested inside a huge target range as complete codistribution;
- it is continuous and avoids selecting a high/low cutoff after outcomes;
- it is computable before any nucleotide identity is opened.

Alternative overlap summaries may be reported descriptively only if frozen before outcome opening; they cannot replace the primary G_st after seeing genetic results.

## 5. Geometry-null centering

Raw T_st can depend on geometry even when source and target spatial processes are biologically unrelated.

Before empirical response opening, simulate lineage-private spatial genetic worlds on the exact fixed dyad geometry using the inherited endpoint-safe IBD and transfer operator.

For each eligible dyad estimate:

- mu0_st = expected T_st under the frozen private-null family;
- optionally sigma0_st = private-null SD where numerically stable.

Primary response is the geometry-centered transfer:

`E_st = T_st - mu0_st`.

This asks whether a source-target pair transfers **more than expected solely from its geometry under independent/private spatial structure**.

A standardized `Z_st` may be a secondary sensitivity if its null SD is stable and prospectively qualified; it cannot replace E_st after empirical opening.

## 6. Primary estimand

Conceptual model:

`E_st = a_s + b_t + beta_G * G_st + beta_P * P_st + beta_X * X_st + error_st`

where:

- a_s absorbs source-specific general usefulness;
- b_t absorbs target-specific predictability;
- beta_G is the primary estimand;
- source-target dependence is handled by a prequalified two-way source/target procedure.

Primary directional hypothesis:

`H1: beta_G > 0`.

Licensed positive interpretation:

> Source-target lineages sharing more of the same landscape show greater excess recurrence of post-IBD spatial genetic structure than expected under lineage-private spatial worlds.

Licensed qualified-null interpretation:

> Within the tested fresh panel and calibrated effect-size domain, geographic co-opportunity did not predict excess source-target recurrence.

A gate failure gives no biological conclusion.

## 7. TTF-Q gate before genetic response

Before any fresh T_st is opened, TTF-Q must characterize G_st on the realized directed-dyad geometry.

Required response-blind outputs:

1. survival after exact source and target identity;
2. additional survival after declared sampling and lineage controls;
3. source/target/dyad signal concentration;
4. null calibration of the intended dyadic inference under the private spatial family;
5. calibration-qualified minimum detectable beta_G across the frozen private-heterogeneity range.

The route stops as `NOT_EVALUABLE_GEOGRAPHIC_RELATION` if the realized G_st is effectively endpoint identity, too concentrated, inferentially miscalibrated, or outside the declared detectable envelope.

## 8. Freshness and data firewall

The new study must exclude every species whose genetic response has already been opened in:

- the completed original fresh 250-species route;
- the conditional two-panel program if any future opening exists;
- the opened host-resource relational study;
- any other exploratory or confirmatory route that exposed pair-specific T_st.

A new deterministic hash namespace must create:

1. development panel: geometry/taxonomy/synthetic worlds only; empirical nucleotide identity remains closed;
2. confirmatory panel: untouched one-shot empirical response after all gates pass.

The source universe remains large enough for this: the prior fresh relational census identified 70,943 metadata candidates after its then-current exclusions, so a genuinely fresh reserve can be selected rather than recycling known responses.

## 9. Secondary lineage-conditioned hypothesis

Only after G_st is frozen as the primary geographic relation should a secondary biological relation R_st be considered.

Preferred structure:

`E_st = ... + beta_G G_st + beta_R R_st + beta_GR (G_st * R_st) + error_st`.

The biologically interesting interaction is beta_GR:

> Does biological similarity determine when shared geography becomes shared phylogeographic structure?

R_st must be chosen from a finite response-blind family before empirical response opening. Candidate domains include dispersal ecology, habitat association and independently reconstructed historical range displacement.

Taxonomic order alone is not privileged; the completed source-pool development showed that order-based conditioning can be geometrically fragile.

## 10. Terminal states

The study has exactly three scientifically valid terminal states:

- `GEOGRAPHIC_COOPPORTUNITY_POSITIVE`: beta_G is positive under qualified inference;
- `GEOGRAPHIC_COOPPORTUNITY_QUALIFIED_NULL`: beta_G does not meet the frozen criterion in a calibrated, detectably informative design;
- `NOT_EVALUABLE_GEOGRAPHIC_RELATION`: response remains closed or uninterpreted because a prerequisite gate fails.

No threshold, overlap metric, subgroup or source-pool rule may be selected after empirical response opening.

## 11. Relation to the place-versus-lineage manuscript

The completed 211-species manuscript remains an unconditional result:

> spatial structure exists within species, but one broad reusable place field was not detected across unseen lineages.

This successor asks the logically next question without changing that result:

> **Is the missing recurrence recovered specifically between lineages that had geographic opportunity to experience the same places?**

A later lineage-relation interaction can then ask:

> **Among codistributed lineages, which biological similarities make the same place produce the same phylogeographic response?**
