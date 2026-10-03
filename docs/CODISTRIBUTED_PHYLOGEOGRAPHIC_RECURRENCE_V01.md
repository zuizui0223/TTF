# Codistributed phylogeographic recurrence — prospective successor v0.1

Status: **response-blind successor; completed 211-species empirical result remains immutable**

## Biological question

The completed genetic TTF result showed that post-IBD spatial genetic structure was detectable within species but that one field learned from 103 species did not predict 108 unseen species across the broad fresh phylogatR panel.

That result does not directly answer the classical comparative-phylogeographic question, because many species in the broad panel do not share the same landscape.

The successor therefore asks a narrower question:

> **When species actually share geographic opportunity, does the spatial configuration of genetic differentiation recur across lineages?**

A second, conditional question asks:

> **Among geographically co-supported species, is recurrence stronger between species belonging to the same broad taxonomic order than between geometry-matched species from different orders?**

The first question is about codistribution. The second tests whether broad lineage similarity conditions the response to the same geography.

## Why this is narrower than a “universal barrier” hypothesis

The response remains post-IBD spatial differentiation, not a labelled physical barrier. A recurrent signal could arise from mountains, rivers, coastlines, climatic transitions, historical fragmentation, refugial structure or other spatially repeated processes.

The successor therefore does not test whether “geographic barriers are universal across animals.” It tests whether **the spatial placement of excess differentiation is reusable among lineages that had an opportunity to encounter the same geography**.

## Separation from the completed 211-species result

The original 211-species result is discovery/motivation only.

No species-level score, nucleotide identity, pairwise genetic distance or favorable subgroup from that result may be used to select the successor panel, thresholds, taxonomic groups or spatial rules.

The successor uses the previously frozen species-disjoint conditional panels:
- 250 development species;
- 250 confirmatory species;
- zero overlap with the original 250-species fresh panel;
- development nucleotide identities remain closed;
- confirmatory nucleotide identities remain closed unless a later formal qualification authorizes opening.

## Geographic co-opportunity

A source species is geographically supported for a target when at least 50% of the target’s graph-edge midpoints lie within 500 km of a source graph-edge midpoint.

At least five source species are required.

This rule was frozen before any successor genetic response was opened.

Response-blind feasibility is adequate:
- development: 109/125 evaluation species have >=5 geographically supported sources;
- confirmatory pre-mask: 105/125;
- confirmatory character-mask survivors: 92/111.

## Geometry-matched lineage contrast

The earlier v0.2 same-order contrast failed Type-I qualification on the confirmatory survivor geometry because same-order and different-order source pools also differed in spatial overlap, centroid distance and graph size. That was an inferential failure, not evidence against the biological hypothesis.

The v0.3 estimator fixes this by matching same-order and different-order sources one-to-one using only response-blind geometry:
- target-to-source coverage;
- source-to-target coverage;
- centroid distance;
- source/target edge-count ratio;
- source/target locality-count ratio.

Every target uses equal numbers of same-order and different-order source species after matching. Source mass is normalized identically on both sides.

Eligible matched targets:
- development geometry: 52;
- confirmatory survivor geometry: 63.

## Development screen

The prospectively frozen v0.3 screen uses four independent private-null reference families (99 worlds each), four private evaluation cells (50 worlds each), and one same-order positive cell (50 worlds).

Inference is a worst-case upper-tail Monte Carlo envelope across the four private reference families.

Frozen screen criteria:
- every private cell: <=5 rejections among 50 worlds;
- same-order A2 positive control: >=35 rejections among 50 worlds.

The development geometry now passes:
- private A0.5: 1/50;
- private A1: 3/50;
- private A2: 2/50;
- private A3: 0/50;
- same-order A2: 40/50.

This is response-blind method development only. It is not a biological result and does not authorize empirical opening.

## Current stopping point

The 214-species confirmatory survivor geometry must pass the same frozen development screen.

The locally recovered survivor geometry reproduces all frozen support counts exactly (92 geographic, 66 same-order, 88 different-order, 63 jointly supported targets), but its byte serialization differs from the canonical frozen survivor CSV. It may be used for numerical development diagnostics but cannot substitute for the canonical byte-bound artifact in a formal authorization.

If both development geometries pass the screen, the next step is to freeze a separately named **formal qualification** with fresh seed namespaces and explicit canonical geometry bindings.

Only a formal qualification PASS may authorize any confirmatory nucleotide identity or genetic transfer response.

## Interpretation branches for a future empirical test

**Positive same-order increment**
: Among species that share geographic support, post-IBD spatial genetic structure is more reusable within broad lineages than across geometry-matched lineages.

**Qualified null**
: Broad taxonomic order does not measurably condition geographic recurrence once co-opportunity and source geometry are controlled.

**NOT_EVALUABLE**
: No biological conclusion.

## Relation to the 211-species manuscript

The flagship manuscript should ask the direct question its frozen estimator answers:

> **Does the spatial configuration of post-IBD genetic differentiation recur in unseen species?**

The broad 211-species result answers that question for an opportunistic multispecies mitochondrial panel. The codistributed successor asks the narrower mechanistic follow-up and should not be folded into the frozen 211-species result unless it completes prospectively.
