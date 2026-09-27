# Butterfly specialization as a hierarchical, multidimensional ecological system

Status: exploratory synthesis of already completed analyses. This document does not introduce a new statistical test.

## Central claim

A butterfly's position on the usual specialist-generalist axis cannot be represented by host-family count alone.

Across the reconstructed S1 ecology panel, five related but non-equivalent ecological dimensions are now distinguishable:

1. **Taxonomic host breadth** — number of larval host-plant families.
2. **Species-level host portfolio richness** — number of resolved host species used within those families.
3. **Geographic resource breadth** — WGSRPD3 geography covered by the union of host distributions.
4. **Anthropogenic resource-opportunity architecture** — how introduced host ranges add geographic opportunity and whether that opportunity is concentrated in one host or distributed across a portfolio.
5. **Climate filtering of realized geography** — how strongly observed butterfly distribution is restricted to climatically closer portions of the reconstructed host-resource envelope.

The main empirical result is not that these dimensions are unrelated. Rather, they are only partially coupled, and host-family breadth compresses ecologically distinct forms of specialization.

## Evidence 1: taxonomic and geographic specialization are not interchangeable

Exact S1 resource reconstruction:

- resource-eligible species: 239
- Spearman(host-family breadth, native geographic resource breadth): **0.276**
- host-taxonomy lower-bound adequate subset: 215 species
- corresponding Spearman: **0.398**
- 33/215 species have absolute taxonomic-vs-geographic rank discordance >= 0.5
- 18 conservative examples are one-family specialists whose geographic resource breadth is in the panel upper quartile

Interpretation:

A butterfly can be taxonomically specialized while possessing a very broad geographic larval-resource base. Family-level specialization therefore does not imply geographic resource specialization.

Source: `benchmarks/exploratory/butterfly_specialization_dimensionality_result_v0.1.json`.

## Evidence 2: human redistribution of hosts expands opportunity across the specialization spectrum

Across 239 resource-eligible butterflies:

- 206/239 (**86.2%**) gain additional reconstructed WGSRPD3 resource opportunity when introduced host ranges are retained
- total native species-units: **26,530**
- total contemporary species-units: **41,083**
- introduced-range addition: **14,553 species-units**
- Spearman(host-family breadth, proportional resource expansion): **0.008**

In the host-taxonomy lower-bound adequate subset:

- 191/215 (**88.8%**) expand
- Spearman(host-family breadth, proportional resource expansion): **0.015**

Median contemporary/native ratios are similar across family-breadth strata:

- 1 family: **1.390**
- 2 families: **1.442**
- 3-5 families: **1.411**
- 6+ families: **1.361**

Interpretation:

Anthropogenic host redistribution does not preferentially increase proportional resource opportunity only for broad generalists. Specialists and generalists can receive comparable proportional expansion.

Source: `benchmarks/exploratory/butterfly_anthropogenic_resource_expansion_result_v0.1.json`.

## Evidence 3: equal aggregate expansion can arise from different host architectures

Among expanded butterflies, total proportional expansion hides a strong structural difference.

Host-taxonomy adequate subset:

- Spearman(host-family breadth, effective contributor number): **+0.492**
- Spearman(host-family breadth, maximum one-host share): **-0.486**
- Spearman(resolved host-species richness, effective contributor number): **+0.776**

Median effective contributor number by host-family stratum:

- 1 family: **1.67**
- 2 families: **2.04**
- 3-5 families: **3.30**
- 6+ families: **5.23**

Median maximum single-host contribution:

- 1 family: **0.747**
- 2 families: **0.566**
- 3-5 families: **0.439**
- 6+ families: **0.313**

Interpretation:

Taxonomic specialists commonly gain anthropogenic opportunity through one or two highly redistributed host species. Broad generalists obtain a similar proportional expansion through a distributed portfolio of many contributing hosts.

Thus host-family breadth predicts **architecture of opportunity** much more clearly than it predicts the total proportional expansion itself.

Source: `benchmarks/exploratory/butterfly_resource_expansion_mechanism_result_v0.1.json`.

## Evidence 4: specialization is hierarchical even within the same family breadth

Among 82 expanded butterflies that all use exactly one host family:

- resolved host-species richness range: **1-37**
- Spearman(resolved host species, effective contributor number): **+0.734**
- Spearman(resolved host species, maximum one-host share): **-0.707**
- Spearman(resolved host species, total introduced-added units): only **+0.269**

One-family host-species richness tertiles:

| host-species richness | median effective contributors | median maximum one-host share |
| --- | ---: | ---: |
| low, median 2 | 1.00 | 1.000 |
| mid, median 4 | 1.74 | 0.692 |
| high, median 12.5 | 2.87 | 0.506 |

The same qualitative relationship occurs within every host-family breadth stratum.

Interpretation:

Two butterflies can both be called "one-family specialists" while having very different ecological portfolios. Family-level specialization and species-level host-portfolio specialization are distinct hierarchical levels.

Source: `benchmarks/exploratory/butterfly_host_specialization_hierarchy_result_v0.1.json`.

## Evidence 5: climate strongly filters realized geography, but not according to host-family breadth

Independent test:

- original independent panel: 32 species
- post-transport-completion pre-climate qualified species: **24**
- climate-informative species: **24**
- median climate-filtering score: **0.801**
- 23/24 species have filtering score > neutral 0.5
- 22/24 species have positive median mismatch difference between never-observed and held-out observed resource units

But the frozen host-breadth release hypothesis was not supported:

- partial Spearman(filtering score, host-family breadth | contemporary resource breadth): **-0.166**
- species-identity-invariant Freedman-Lane residual permutations: **9,999**
- one-sided p = **0.2237**
- decision: **PILOT_DERIVED_CLIMATE_RELEASE_HYPOTHESIS_NOT_SUPPORTED**

Median filtering scores overlap strongly across host-family strata:

- 1 family: **0.796**
- 2 families: **0.748**
- 3-5 families: **0.813**
- 6+ families: **0.801**

Interpretation:

Climatic filtering within available larval-resource geography is common, but broad taxonomic diet breadth does not detectably release butterflies from that filtering.

Source after climate PR integration:
`benchmarks/exploratory/butterfly_climate_release_postgate_independent_result_v0.1.json`.

## Integrated ecological interpretation

The usual specialist-generalist label mixes several mechanisms that should be separated.

A useful hierarchy is:

```
host-family breadth
        |
        v
species-level host portfolio
        |
        +------> geographic resource breadth
        |
        +------> architecture of anthropogenic opportunity
                         |
                         v
              potential resource geography
                         |
                         v
             climate + dispersal + life history
                         |
                         v
                realized butterfly geography
```

The data support partial coupling between these levels, not a single latent specialization axis.

In particular:

- Family breadth is a coarse taxonomic descriptor.
- Host-species richness better resolves portfolio structure within a family-breadth category.
- Geographic resource breadth can be broad even in family specialists.
- Human redistribution of hosts expands opportunity broadly across specialists and generalists.
- Specialists and generalists reach similar aggregate expansion through different portfolio architectures.
- Climate filtering is strong across many butterflies but is largely orthogonal to host-family breadth.

## Paper-level thesis

**Butterfly specialization is hierarchical and multidimensional: family-level diet breadth masks species-level host portfolios that structure anthropogenic resource opportunity, while climate independently filters the realized portion of that opportunity.**

This framing converts the earlier negative results into ecological information rather than treating them as failed predictors.

The important result is not merely that host-family breadth sometimes performs poorly. It is that each apparent failure identifies a different ecological layer that the one-dimensional specialist-generalist label collapses.

## Suggested figure sequence

1. **Taxonomic vs geographic specialization**  
   Host-family breadth against native resource geography; highlight specialist-wide discordant species.

2. **Anthropogenic expansion without generalist advantage**  
   Native vs contemporary resource opportunity and proportional expansion across host-family strata.

3. **Different architectures produce similar expansion**  
   Effective contributor number and maximum one-host share across host-family breadth.

4. **Hierarchy within specialists**  
   One-family species: resolved host-species richness against effective contributor number / single-host dominance.

5. **Climate as a separate filter**  
   Independent climate-filtering score across family-breadth strata, with the frozen non-supported host-breadth test.

## Claim boundaries

This synthesis remains ecological opportunity analysis, not a claim that butterflies actually use every introduced host population.

Specifically:

- WCVP introduced status is a geographic distribution status, not an introduction date.
- HOSTS completeness is not guaranteed.
- WGSRPD3 is coarse geography.
- A reconstructed host-resource unit is potential opportunity, not realized larval use.
- A never-observed butterfly unit is not a confirmed true absence.
- Climate filtering is inferred using the frozen cross-fit and sampling-effort proxy, not occupancy surveys.
- The species-level portfolio-to-climate connection remains a post-result hypothesis and should require a fresh panel before confirmatory testing.
