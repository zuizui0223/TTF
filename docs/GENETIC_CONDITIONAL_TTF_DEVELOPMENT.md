# Conditional genetic TTF successor study

Status: **response-blind successor frozen before any nucleotide identity from the new panels is opened**.

The completed 211-species phylogatR result remains an immutable unconditional baseline. Its empirical conclusion is not being repaired or re-tested. The successor asks a different question:

> **Among species that share geographic support, does broad biological similarity make spatial genetic structure more transferable?**

The motivation is straightforward: arbitrary species differ in movement, demography, historical range dynamics, habitat use and barrier response. At the same time, the original TTF field already uses a 500-km Gaussian kernel, so distant training edges are strongly down-weighted automatically. Geographic mismatch therefore cannot simply be assumed to explain the completed non-transfer result.

The successor consequently uses geographic co-coverage as a **support control**, then tests coarse lineage similarity as the primary increment.

## Frozen two-panel design

The exact authenticated phylogatR archive is reused only as a source universe. The original 250-species fresh panel is excluded completely. Among the remaining response-blind eligible species, one fixed SHA-256 ranking supplies two non-overlapping panels:

- **development panel:** ranks 1–250; geometry/taxonomy and synthetic development only; its nucleotide identities remain closed permanently;
- **confirmatory panel:** ranks 251–500; reserved for one future empirical opening after qualification.

The response-blind census found 14,612 eligible species after exclusions and no geometry failures among the first 500 selected species. The development and confirmatory panels have zero species overlap with each other and with the earlier fresh panel.

## Geographic support control

For each held-out target species, a training species is geographically supported only when at least **50% of the target's graph-edge midpoints lie within 500 km of a source graph-edge midpoint**. At least five source species are required.

This support rule is not the main biological hypothesis. It creates a matched geographic opportunity set before lineage similarity is tested.

To prevent a smaller source pool from changing the comparison merely by reducing total training weight, every target-specific field rescales selected-source edge weights by `N_train / N_selected`. Thus the unconditional, geography-conditioned and geography+order fields retain the same total pre-kernel training-species weight mass; only source composition and its realized spatial kernel geometry change. The synthetic qualification calibrates the remaining finite-geometry behavior.

Feasibility is adequate before outcomes: 109/125 development targets and 105/125 confirmatory targets have at least five geographically supported source species.

A paired diagnostic remains available:

[
\Delta_{geo,t}=C_{geo,t}-C_{all,t}.
]

Because the inherited TTF field is already Gaussian-local at 500 km, this diagnostic asks only whether explicit source-support restriction adds predictive value beyond locality weighting already built into TTF.

## Primary hypothesis — same-order transfer beyond geography

The primary field starts from the same geographic support rule and additionally requires source and target to share the exact non-empty phylogatR taxonomic **order**.

For each supported target,

[
\Delta_{order,t}=C_{geo+order,t}-C_{geo,t}.
]

The inferential statistic is the species-equal mean of this paired increment.

This is deliberately stronger than asking whether one favorable order has a positive TTF score. The same rule is applied to every response-blind supported target, and the comparison is within target against a geography-matched reference field.

Seventy-two of 125 development targets and 78 of 125 confirmatory targets already have at least five geographically supported same-order source species.

Family-level conditioning is excluded prospectively because response-blind support is too sparse.

## What taxonomic order does and does not mean

Order is a coarse pre-outcome proxy for shared biology. It may correlate with dispersal mode, body plan, life history or historical response, but this study does **not** claim that order measures any of those mechanisms directly.

A positive primary result would mean:

> among species with comparable geographic support, spatial genetic patterns are more reusable when source and target belong to the same broad lineage.

It would not identify why.

Direct dispersal traits, generation time, reconstructed range history or barrier exposure require a separately frozen external-data layer.

## Inference and qualification

The primary same-order increment and the geographic diagnostic use one-sided centered, studentized bootstrap inference over held-out species (1,999 resamples; alpha 0.05). Species remain the inferential units.

Before any confirmatory nucleotide identity can be opened, the estimator must pass synthetic Type-I and power qualification on the development geometry. The confirmatory panel must then pass character-mask admissibility and exact survivor-geometry requalification.

The development panel's empirical nucleotide identities remain unopened permanently.

## Interpretation

The key future branches are:

- **same-order increment positive:** lineage similarity adds transferable information beyond geographic co-coverage;
- **same-order increment null:** broad taxonomic similarity does not measurably rescue transfer once geography is matched;
- **geographic diagnostic positive:** explicit range-support restriction adds value beyond the already-local kernel;
- **geographic diagnostic null:** the inherited Gaussian field already handles geographic mismatch adequately.

None of these outcomes identifies dispersal ability or historical range dynamics by itself.
