# Co-distribution-conditioned place recurrence v0.1

Status: **new response-blind successor study; frozen before any genetic response from the 1,000-species universe is opened**

## Question

The completed 211-species result asked whether one broad post-IBD geographic field transfers across arbitrary unseen lineages. It did not. That result leaves a narrower comparative-phylogeographic question:

> **If two lineages actually occupy and are sampled across the same geography, does the location of post-IBD genetic differentiation recur across them?**

This is deliberately narrower than asking whether geographic barriers are universal across animals.

## Why this is a new study rather than a rescue

The original 211-species empirical result is immutable. The conditional-order v0.2/v0.3 route is also closed and keeps inferential alpha zero. Relational v0.3 is closed. No old test is reopened.

This successor uses a 1,000-species genetic-response universe that excludes all 250 original fresh species and all 500 species from the earlier conditional two-panel design. Its geography and taxonomy have been inspected response-blind, but nucleotide identity, pairwise genetic distances and genetic transfer outcomes remain unopened.

The new test does **not** compare one source pool against another. It first defines a biologically coherent analysis population from geography alone, then asks whether place recurrence exists within that population.

## Co-distribution definition

For species A and B, let C(A->B) be the fraction of A's frozen graph-edge midpoints whose nearest B edge midpoint lies within 500 km.

Define symmetric co-distribution as:

`C_sym(A,B) = min(C(A->B), C(B->A))`.

The primary rule is `C_sym >= 0.50`. Thus at least half of the sampled edge geometry of **both** species must be geographically represented by the other species. A held-out target requires at least five such training species.

The 500-km radius is inherited from the established TTF spatial field; the 0.50 threshold and five-source minimum are frozen before genetic response opening.

## Response-blind feasibility

The source universe contains 1,000 species and 167,828 frozen edge midpoints. Across all 1,000 species, 56,818 undirected pairs satisfy the primary symmetric co-distribution rule.

A new hash namespace splits the species into 500 development and 500 confirmatory species, then into 250 training and 250 evaluation species inside each panel.

**Development:** 200/250 evaluation species have at least five strongly codistributed training sources; the median is 33 sources. Supported targets span 17 orders, with Lepidoptera 94/200 (47%).

**Confirmatory:** 203/250 evaluation species have at least five strongly codistributed training sources; the median is 32.5 sources. Supported targets span 23 orders, with Lepidoptera 90/203 (44.3%).

So the question is empirically supportable without turning into a Lepidoptera-only or same-order analysis.

## Primary estimand

For each supported held-out species:

1. remove ordinary within-species IBD with the inherited endpoint-safe leave-two-localities-out procedure;
2. learn a 500-km Gaussian geographic field only from its strongly codistributed training species;
3. predict the target's post-IBD edge turnover;
4. score prediction with Spearman correlation;
5. average scores equally across held-out species.

Selected-source weights are scaled by `N_train/N_selected`, as in the prior target-conditioned implementation, so the total pre-kernel training mass and prior influence do not change merely because a target has fewer eligible source species.

There is **no primary comparison to all-species, same-order, different-order, host-similar or climate-similar fields**.

## Hypotheses

**Place-recurrence hypothesis:** once shared geographic opportunity is required prospectively, the same places contain reusable information about post-IBD differentiation in unseen lineages.

**Lineage-conditioned hypothesis:** even among strongly codistributed lineages, spatial genetic structure remains predominantly lineage-conditioned, so the held-out geographic recurrence statistic is not extreme relative to private spatial structure.

## Qualification before any empirical response

Development uses synthetic genetic worlds only. The estimator must control false recurrence under private lineage-specific spatial fields across residual amplitudes A=0.5, 1, 2 and 3, and recover a frozen shared A=2 positive field.

The confirmatory panel then receives character-mask admissibility only. After mask failures are removed without replacement or graph repair, the exact survivor geometry is requalified. Only a passing survivor geometry can authorize nucleotide-identity opening.

The three terminal states are:

- **PLACE_RECURRENCE_DETECTED**
- **QUALIFIED_NONDETECTION_WITHIN_CODISTRIBUTED_LINEAGES**
- **NOT_EVALUABLE**

## Biological interpretation

A positive result would show that comparative phylogeographic concordance becomes predictively reusable once species are required to share the same sampled geography.

A qualified non-detection would be stronger than the original 211-species result: it would show that geographic mismatch is not sufficient to explain the lack of recurrence, because recurrence still fails among lineages with substantial mutual geographic coverage.

Neither branch identifies a named river, mountain, refugium or ecotone as causal. The response is post-IBD spatial differentiation, not a direct barrier label.

## Relation to comparative phylogeography

This design targets the classic **codistributed-species** problem rather than an animal-wide universal barrier map. Its novelty is the prospective prediction criterion: a place component counts as recurrent only if geography learned without the target species predicts that unseen target.

The study therefore separates two questions that are often conflated:

1. do observed species show apparently concordant spatial breaks?
2. is that concordance strong enough to predict an unseen codistributed lineage?

Only the second is tested here.

## Hard boundaries

No subgroup genetic-response search is authorized. Order, family, habitat and climate can be reported descriptively for panel composition, but cannot be selected as rescue moderators after the primary result. Any later mechanistic analysis requires a separately named fresh hypothesis and independent response allocation.
