# Genetic TTF — place-versus-lineage reframing v0.1

Status: **alternative ecological framing on a new branch; frozen estimands, thresholds, split and empirical results unchanged**

Branch: `paper/place-vs-lineage-phylogeography-v0-1`

## Biological problem

Comparative phylogeography has traditionally used concordant genetic discontinuities among species as evidence that geography, landscape barriers or shared history structure genetic variation repeatedly across lineages.

The stronger question is:

> **When a species shows excess genetic differentiation at a particular place beyond ordinary isolation by distance, does that geographic configuration recur strongly enough to predict differentiation in an entirely unseen lineage?**

This reframes the frozen TTF estimand as a test of **geographic recurrence across lineages**, not as an abstract transfer problem.

## Core contrast

### Place-recurrence hypothesis

If landscape position has a lineage-general component — for example through mountain systems, rivers, coastlines, climatic transition zones, historical refugial boundaries, or other persistent spatial constraints — then post-IBD differentiation learned from training species should predict where unseen species show relatively strong differentiation.

Operational prediction:

[
T = operatorname{mean}_s ho(widehat{D}^{,postIBD}_s, D^{postIBD}_s) > 	ext{qualified null expectation}.
]

A positive result would establish a reusable place component. It would not, by itself, identify a specific barrier mechanism.

### Lineage-conditioned hypothesis

If the effect of geography depends strongly on dispersal, habitat association, demographic history, range dynamics, barrier permeability or other lineage-specific properties, then within-species spatial structure may remain detectable while its geographic configuration fails to recur across unseen species.

Operational prediction:

- within-species self structure is detectable under its separately qualified reference;
- cross-species post-IBD recurrence is not detected.

This is the branch observed in the frozen empirical result.

## Frozen empirical evidence

Fresh phylogatR COI/COX1 panel:

- 211 Phase-2 survivor species;
- 103 training species;
- 108 entirely unseen evaluation species;
- exact-geometry shared-signal qualification passed before nucleotide identity opening;
- shared-A2 power lower bound: 0.971387;
- post-IBD cross-species statistic: (T=0.0325655), profiled-private (p=0.739261);
- separately qualified within-species self test: (p=0.008991);
- frozen decision: `LINEAGE_CONDITIONED_SPATIAL_STRUCTURE_WITHIN_TESTED_DOMAIN`.

The key empirical asymmetry is therefore:

> **spatial genetic structure is detectable within species, but its geographic configuration is not detectably recurrent across unseen lineages under the tested spatial field.**

## Why this is a comparative-phylogeography question

The framing connects directly to the field's long-standing concordance problem.

Traditional logic:

[
	ext{similar breaks in observed taxa}
Rightarrow
	ext{shared biogeographic process}.
]

Predictive logic used here:

[
	ext{shared geographic component}
Rightarrow
	ext{structure learned without species }s
	ext{ should predict species }s.
]

The second criterion is intentionally stronger because the evaluation lineage is absent from field fitting.

This complements trait-based comparative phylogeography. Trait-based studies ask why some taxa develop stronger or weaker phylogeographic structure. This analysis asks whether the **location** of post-IBD structure itself is reusable across lineages.

## Geographic barriers: correct claim boundary

Geographic barriers are a major biological motivation, but **the response is not a direct barrier annotation**.

Post-IBD spatial differentiation may reflect:

- physical barriers;
- ecological transition zones;
- historical range fragmentation;
- refugial structure;
- secondary contact;
- spatially heterogeneous demography;
- combinations of these processes.

Therefore the manuscript may say:

> If major geographic or historical constraints impose similar spatial effects across many lineages, they should contribute to a recurrent geographic component that is predictable in unseen species.

It must not say:

> The analysis proves that geographic barriers are absent, unimportant, or non-general across animals.

## Scope of the 211-species panel

The panel is a broad, opportunistically assembled phylogatR COI/COX1 animal dataset after frozen geometry and character-support filters. It is not a probability sample of Animalia. Insects and other arthropods are prominent, with smaller representation from other animal groups.

Preferred scope language:

- "a broad multispecies animal COI/COX1 panel";
- "211 animal species represented in the fresh phylogatR confirmatory domain";
- "across the tested multispecies mitochondrial panel".

Avoid:

- "all animals";
- "animal-wide universal rule";
- "global proof that barriers are lineage specific".

## Main claim ladder

### Directly supported

1. The exact 211-species geometry was qualified to distinguish private spatial structure from the declared shared positive field before empirical genetic outcomes were opened.
2. Within-species spatial structure was detectable relative to the separately frozen structural-null reference.
3. A reusable post-IBD geographic field learned from 103 species did not predict differentiation in 108 unseen species under the qualified test.
4. The observed combination is consistent with lineage-conditioned spatial structure within the tested domain.

### Mechanistically plausible but not identified

- lineage differences in dispersal;
- habitat association;
- demographic history;
- historical range movement;
- permeability of the same physical feature;
- taxon-specific response to climatic or geological history.

### Not supported by the current test

- absence of geographic barriers;
- proof that every barrier is lineage specific;
- identification of any particular mountain, river, refugium or transition zone as causal;
- whole-genome generalization;
- universality across taxa or spatial scales.

## Recommended manuscript question

> **Is phylogeographic structure a property of place or lineage?**

Operational subtitle/question:

> **Does post-IBD spatial differentiation recur strongly enough across species that geography learned from some lineages predicts unseen lineages?**

## Candidate titles

1. **Is phylogeographic structure a property of place or lineage? A predictive test across 211 mitochondrial datasets**
2. **Does phylogeographic structure recur across species? A prospective test of geographic recurrence in 211 animal datasets**
3. **Spatial genetic structure is detectable within species but not predictably recurrent across lineages**
4. **From concordance to prediction: testing the geographic recurrence of phylogeographic structure across 211 species**

Preferred working title: **Is phylogeographic structure a property of place or lineage? A predictive test across 211 mitochondrial datasets**

## Novelty statement

The novelty is not the observation that species can respond differently to shared barriers. Comparative phylogeography has long recognized concordance and discordance.

The new contribution is to turn geographic recurrence into an **out-of-lineage prediction criterion**:

> a geographic component is considered reusable only if it predicts the spatial configuration of differentiation in species excluded entirely from model fitting, after ordinary IBD is removed and after the exact sampling geometry has been qualified for false-transfer control and shared-signal detectability.

## Relation to TTF methodology

TTF remains the inferential machinery, not the biological headline.

The ecological narrative is:

[
	ext{shared barriers/history}
ightarrow
	ext{recurrent place-specific differentiation}
ightarrow
	ext{prediction in unseen lineages}.
]

TTF provides the prospective split, endpoint-safe IBD residualization, exact-geometry qualification and abstention rules needed to test that narrative without retrospective concordance.

## Next-step research question

The frozen result motivates, but does not answer:

> **Which pairs of lineages should share a geographic barrier response?**

A future fresh study could prospectively test whether recurrence is conditional on dispersal traits, habitat, present climatic niche or shared historical range displacement. That question should be treated as a new hypothesis family rather than used to retune the completed 211-species result.
