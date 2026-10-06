# Codistributed recurrence — literature positioning v0.1

Status: **post-result positioning; frozen estimand and interpretation branch preserved**

## The field-level problem

Comparative phylogeography has historically treated concordance among codistributed taxa as evidence that shared landscape history, geographic barriers, or climatic events repeatedly structured genetic variation. Classic analyses therefore compare already-observed phylogeographic trees, breaks, or regional partitions across taxa.

That literature establishes the biological importance of the question, but it does not make the same inferential claim as the present prospective design.

## Closest precedents

### Lapointe & Rissler (2005)

Lapointe, F.-J. & Rissler, L.J. (2005). Congruence, consensus, and the comparative phylogeography of codistributed species in California. *The American Naturalist* 166:290–299. DOI: 10.1086/431283.

They tested concordance among nine already-observed phylogeographic datasets and synthesized them into a regional supertree. This is a direct precedent for asking whether codistributed species share geographic signal.

**Difference here:** the present study asks whether spatial information learned from source lineages predicts the post-IBD differentiation pattern of a target lineage that is excluded from fitting.

### Papadopoulou & Knowles (2016)

Papadopoulou, A. & Knowles, L.L. (2016). Toward a paradigm shift in comparative phylogeography driven by trait-based hypotheses. *PNAS* 113:8018–8024. DOI: 10.1073/pnas.1601069113.

They argued that a strict focus on concordance can make tests too generic and can lead to ad hoc explanations of discordance. Their proposed shift is toward explicit predictions about why species respond differently to shared geography.

**Difference here:** the present study first asks a more primitive question — whether the amount of genuinely shared sampled geography itself predicts how reusable a lineage's spatial genetic pattern is, before invoking a particular trait mechanism.

### Massatti & Knowles (2014)

Massatti, R. & Knowles, L.L. (2014). Microhabitat differences impact phylogeographic concordance of codistributed species: genomic evidence in montane sedges from the Rocky Mountains. *Evolution* 68:2833–2849. DOI: 10.1111/evo.12491.

This work demonstrates that broadly codistributed taxa can respond differently to the same glacial landscape because microhabitat changes effective historical connectivity.

**Relevance here:** it supplies a concrete mechanism for lineage-conditioned geography. The same broad place need not have the same effective permeability for every lineage.

### Edwards et al. (2022)

Edwards, S.V. et al. (2022). The evolution of comparative phylogeography: putting the geography (and more) into comparative population genomics. *Genome Biology and Evolution* 14:evab176. DOI: 10.1093/gbe/evab176.

This review describes comparative phylogeography as explicitly place-based: landscape features such as mountains, rivers and transition zones are sought as drivers of genomic splits shared among species.

**Relevance here:** the completed 211-species TTF result and the new codistribution study turn that shared-place idea into an out-of-lineage prediction criterion.

### Jensen et al. (2024)

Jensen, A.J. et al. (2024). Geographic barriers but not life history traits shape the phylogeography of North American mammals. *Global Ecology and Biogeography* 33:e13875. DOI: 10.1111/geb.13875.

Jensen et al. synthesized 229 North American mammals, mapped published genetic breaks, and tested whether known break locations were associated with mountains, major water bodies, ecoregion boundaries and terrain. They found barrier associations but substantial among-species heterogeneity, including within orders.

This is the closest large-scale empirical precedent.

**Critical difference:** Jensen et al. condition on break locations already inferred for each focal species and then model which geographic features align with those breaks. The present study does not require a target's break map during field construction. Instead it asks whether a source lineage's post-IBD spatial information is reusable for an entirely unseen target, and whether that reuse increases continuously with response-blind geographic co-opportunity.

## Novel inferential contrast

The classical concordance question can be written schematically as:

`Given observed phylogeographies of species 1...n, how concordant are their breaks or histories?`

The prospective recurrence question is:

`Given source lineages only, does their spatial information predict an excluded target lineage, and is that prediction stronger when source and target genuinely share more sampled landscape?`

Those are not equivalent.

Retrospective concordance can be scientifically valuable even when it would not generalize to an unseen lineage. Conversely, successful out-of-lineage prediction is direct evidence that a component of the spatial pattern is reusable rather than merely recognizable after all outcomes are visible.

## Why G_st is not simply "distance"

The frozen focal relation is

`G_st = min(target->source 500-km edge-midpoint coverage, source->target 500-km edge-midpoint coverage)`.

It represents mutual sampled-landscape opportunity, not centroid proximity.

The response-blind TTF-Q audit already shows that G_st retains substantial independent information after exact source identity, exact target identity, edge-count imbalance, locality-count imbalance, same-order membership and centroid distance:

- development final unique fraction: 0.2149;
- confirmatory final unique fraction: 0.2054;
- maximum source or target signal concentration <2%.

Thus the prospective predictor is not reducible to the statement that nearby species have nearby centroids.

## The mechanical-opportunity problem and its solution

Greater geographic overlap can mechanically make a transfer operator better supported. A naive positive association between G_st and raw T_st would therefore be ambiguous.

The frozen v0.3 fixed-dyad design addresses this explicitly:

1. run the **same** endpoint-safe post-IBD single-source TTF operator in lineage-private synthetic worlds;
2. estimate each dyad's geometry-specific expected transfer, `mu0_st`;
3. analyze `E_st = T_st - mu0_st`;
4. calibrate beta_G against independent private-world reference distributions;
5. take the least-favourable Monte Carlo p-value across the frozen private-amplitude family.

Therefore a confirmatory positive result requires more than "overlap makes the kernel estimable." It requires a G_st slope exceeding what the same sampling geometry produces when every lineage has private spatial structure.

## Observed result in the comparative-phylogeographic context

The one-shot confirmatory test followed the prospectively frozen qualified-null branch.

- empirical `beta_G = -0.004447`;
- least-favourable Monte Carlo envelope `p = 0.969`;
- survivor synthetic qualification had already shown calibrated private-null Type-I control and **0.946** power for the frozen shared-place positive control.

Thus the completed broad non-transfer result is not rescued by explicitly conditioning on how much sampled geography two insect lineages share. The important point is not merely another failure to find concordance. The study prospectively tested the most direct geographic-opportunity explanation for discordance and rejected it under a design that had sufficient response-blind relation information and calibrated power.

This sharpens the classical comparative-phylogeographic problem. Shared geography can matter strongly within individual lineages while still failing to define a portable spatial genetic field across lineages. The same place can therefore have different effective phylogeographic meaning for different organisms.

The next mechanistic level is not a search for a single hidden universal barrier map. It is the interaction between place and lineage: dispersal, habitat dependence, historical occupancy, demographic response, and other processes that change how a geographic feature is translated into genetic structure.

## Claim hierarchy

### A positive result would support

> Among the tested fresh insect lineages, greater mutual geographic co-opportunity is associated with greater **excess cross-lineage recurrence** of post-IBD spatial genetic structure than expected from lineage-private spatial worlds.

This would provide a predictive bridge between classical codistribution and phylogeographic concordance.

### The observed qualified null supports

> Even among lineages with varying and often substantial shared geographic opportunity, greater co-distribution did not generate detectable excess recurrence at the tested scale.

Combined with the completed 211-species result, this shows that simple failure to occupy or sample the same geography is insufficient to explain the weak cross-lineage transfer in the qualified insect domain.

### Neither result directly establishes

- that a named mountain, river, coastline or refugium is causal;
- that all barriers are lineage-specific;
- that mitochondrial structure generalizes to nuclear genomic structure;
- that the result is universal across Animalia.

## Recommended manuscript positioning

The paper should not sell the contribution as "the first study of concordance among codistributed species." That would be false.

The defensible novelty is:

> **a prospective, out-of-lineage test showing that greater shared geographic opportunity is not sufficient to make phylogeographic structure reusable across lineages, after calibrating the measurement opportunity itself with lineage-private spatial null worlds.**

This moves the comparative-phylogeographic problem from asking only whether breaks are concordant to asking why the same landscape is translated into different genetic structure by different lineages.
