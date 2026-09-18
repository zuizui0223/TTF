# Conditional genetic TTF successor study

Status: **response-blind successor frozen before any nucleotide identity from the new panels is opened**.

The completed 211-species phylogatR result remains an immutable unconditional baseline. Its empirical conclusion is not being repaired or re-tested. The successor asks a different question:

> **When does a spatial genetic pattern transfer across species?**

The biological motivation is that absolute geographic concordance should not be expected across arbitrary species with different distributions, dispersal, demographic histories and barrier responses. A useful cross-species test therefore needs to distinguish lack of shared geographic opportunity from lack of transfer conditional on shared opportunity.

## Frozen design

The exact authenticated phylogatR archive is reused only as a source universe. The original 250-species fresh panel is excluded completely. Among the remaining response-blind eligible species, a fixed SHA-256 ranking supplies two new non-overlapping panels:

- **development panel:** ranks 1–250; geometry/taxonomy and synthetic development only; its nucleotide identities remain closed permanently;
- **confirmatory panel:** ranks 251–500; reserved for one future empirical opening after qualification.

The response-blind census found 14,612 eligible species after exclusions and no geometry failures among the first 500 selected species. The development and confirmatory panels have zero species overlap with each other and with the earlier 250-species fresh panel.

## H1 — geographic conditioning

For each held-out target species, a training species enters the conditional field only when at least **25% of the target's graph-edge midpoints lie within 500 km of a source graph-edge midpoint**. At least five source species are required.

The primary statistic is paired within target:

[
\Delta_{geo,t} = C_{geo,t} - C_{all,t},
]

where (C_{all,t}) is the unchanged unconditional TTF score and (C_{geo,t}) uses only the prospectively defined co-covered source pool. The panel statistic is the species-equal mean of (Delta_{geo,t}) over response-blind supported targets.

This directly tests whether the earlier unconditional analysis was asking too much by pooling species that do not share the relevant geographic support.

Feasibility is adequate before outcomes: 113/125 development targets and 108/125 confirmatory targets have at least five qualifying source species.

## H2 — lineage similarity beyond geography

The secondary field applies the **same geographic support rule** and additionally requires source and target to share the exact non-empty phylogatR taxonomic **order**. The statistic is

[
\Delta_{order,t} = C_{geo+order,t} - C_{geo,t}.
]

Thus order is tested only as an increment beyond geographic co-coverage, not by cherry-picking an order after viewing genetic responses. Seventy-nine of 125 development targets and 86 of 125 confirmatory targets already have at least five geographically supported same-order source species.

Family-level conditioning is excluded prospectively because response-blind support is too sparse. Taxonomic order is only a coarse proxy for shared biology; it is not a direct measure of dispersal ability, generation time or historical range dynamics.

## Inference

Both increments use a one-sided centered, studentized bootstrap over held-out species (1,999 resamples; alpha 0.05). Species remain the inferential units.

Before any confirmatory nucleotide identity can be opened, the new estimators must pass separate synthetic Type-I and power qualification on the development geometry, followed by exact-geometry requalification after confirmatory character-mask admissibility.

## Interpretation

Possible outcomes are now biologically sharper than the unconditional result:

- **H1 positive:** spatial genetic transfer exists when species genuinely share geographic support;
- **H1 null:** co-distribution alone does not rescue cross-species transfer;
- **H2 positive:** coarse lineage similarity adds predictive value beyond co-distribution;
- **H2 null:** same-order membership does not add detectable transfer once geography is matched.

None of these outcomes by itself identifies dispersal ability or distributional history as the mechanism. Direct traits or reconstructed histories would require a new prospectively frozen data-integration layer.
