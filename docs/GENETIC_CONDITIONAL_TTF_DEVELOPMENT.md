# Conditional genetic TTF development

Status: **response-blind development on a species-disjoint unused phylogatR panel; nucleotide identity unopened**.

## Why this successor exists

The completed genetic TTF one-shot result is an unconditional global baseline. It showed that a geographic field learned from a broad set of species did not predict post-IBD differentiation in unseen species under the frozen global field. That result remains immutable.

The biological limitation is straightforward: different species differ in range history, dispersal, habitat and barrier response. A universal absolute map is therefore a deliberately strong null model of sharedness. The successor question is not whether the old result can be rescued. It is:

> **When does phylogeographic structure become transferable across species?**

The first condition must be geographic: species cannot meaningfully share an absolute spatial pattern where their sampled ranges do not overlap. Only after that support requirement is fixed do we test whether lineage similarity improves transfer.

## New species-disjoint panel

The exact authenticated phylogatR archive is reused as a source database, but every species in the completed parent Phase-1 panel is excluded before selection.

Response-blind census:

- exact source SHA-256: `5a0fd9ac25893c749d14186fbcce4a46b99163c9d810b36e40eebce7bece61a5`;
- prior Phase-1 species excluded: **250**;
- eligible species after the prior-panel exclusion: **14,612**;
- new selected panel: **250 species**;
- overlap with the old 250: **0**;
- first 250 new hash ranks passing geometry: **250/250**;
- deterministic split: **125 training / 125 evaluation**;
- nucleotide identity opened: **no**.

This is a species-disjoint prospective panel within the same source database, not an independent-database replication.

## Geographic support definition

For source species `s` and target species `t`, directed co-distribution is

> the fraction of target edge midpoints lying within **500 km** of at least one source edge midpoint.

A source is called geographically supported for the target when this fraction is at least **0.25**.

The 500-km radius is inherited from the frozen genetic TTF field scale. The 0.25 threshold and minimum of three sources were fixed using geometry/taxonomy only, before any conditional-panel nucleotide identity or genetic outcome was opened.

Under the new 125/125 split:

- targets with at least 3 geographically supported sources: **113 / 125**;
- median number of geographically supported sources: **33** (IQR 20–39).

## Primary biological test: lineage similarity given geography

Taxonomic order is used only as a coarse, response-blind **lineage-similarity proxy**. It is not labelled dispersal ability and it is not treated as evidence of common demographic history.

For targets having at least three geographically supported same-order training species:

- supported evaluation targets: **80 / 125**;
- median same-order geographically supported sources: **6** (IQR 2–13).

For each fixed target:

1. **Geo field:** learn from every geographically supported training species.
2. **Geo + order field:** learn only from the geographically supported training species in the same taxonomic order.
3. Score both fields on the same target species.
4. Compute the paired target-level increment.

The primary estimand is

[
\Delta_{O|G}
=
\frac{1}{|S_*|}
\sum_{t\in S_*}
\left(
C_{t,\,\mathrm{same\ order+geo}}
-
C_{t,\,\mathrm{geo}}
\right).
]

The scientific hypothesis is `Delta_order_given_geo > 0`.

A positive qualified result would mean that lineage similarity improves spatial transferability **among species that actually share geographic support**. It would not identify dispersal, life history, barrier permeability or shared history as the mechanism.

## Secondary geometry diagnostic

A secondary development contrast is

[
\Delta_G
=
\mathrm{mean}
\left(
C_{\mathrm{geo}} - C_{\mathrm{all}}
\right)
]

on the 113 targets with at least three geographically supported sources.

This asks whether excluding globally irrelevant species helps prediction. Because the TTF kernel is already spatially local, this contrast is treated as a development diagnostic until a dedicated geography-conditioned synthetic positive-control family is frozen and qualified.

## Why family is not tested

Family-level similarity was evaluated only for response-blind support feasibility.

Only **11 / 125** evaluation targets have at least three family-matched training species meeting the 25% geographic-support rule; **87 / 125** have none. Family is therefore rejected as a confirmatory condition **before genetic outcomes** rather than searched post hoc.

## Synthetic qualification before biology

No conditional empirical opening is authorized yet.

The primary paired increment must first pass exact-geometry synthetic qualification:

- **private null:** every species receives its own independent spatial boundary field;
- **same-order positive control:** species in the same frozen taxonomic order share a boundary field, while different orders receive independent fields;
- positive-control residual amplitude: **2.0**;
- noise SD: **0.10**;
- IBD strength: **1.0**;
- alpha: **0.05**;
- private-null Wilson 95% upper rejection bound: **<= 0.10**;
- positive-control Wilson 95% lower power bound: **>= 0.80**.

The candidate empirical inference is a one-sided centered held-out-species bootstrap of the **paired per-target increments**. It is licensed only if the synthetic gate shows that this inference controls Type I error and has adequate power on the exact new geometry.

Failure means `NOT_EVALUABLE` and conditional nucleotide identity remains closed.

## Dispersal and life-history extension

The user's biological objection is broader than taxonomy: movement ability and distribution history differ among species.

Those variables are scientifically important, but they are **not present as direct measurements in the frozen phylogatR panel**. They may enter only through a separately frozen external trait source with:

1. an explicit trait definition;
2. a response-blind species-matching rule;
3. a predeclared missing-data rule;
4. adequate coverage on the new panel;
5. no genetic outcome used to choose traits or thresholds.

Until then, order is described only as lineage similarity, never as measured dispersal ability.

## Current state

The response-blind panel and support census are feasible and frozen. Generic conditional source restriction, paired increments and grouped synthetic genetic worlds are implemented. No conditional-panel nucleotide identity, pairwise genetic distance or empirical conditional TTF result has been opened.
