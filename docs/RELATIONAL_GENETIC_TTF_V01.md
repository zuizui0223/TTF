# Relational genetic TTF v0.1 — fresh hypothesis contract

Status: **prospective successor; response must remain closed until the relational design is frozen and qualified**.

## Separation from completed TTF studies

This study does not rescue, retune, reinterpret, or rerun the completed unconditional phylogatR result. It also does not modify the conditional same-order v0.2 design after its survivor-geometry Type-I failure. Those results remain immutable.

The new question is:

> **Is cross-species transferability of spatial genetic structure predictable from the ecological relationship between a source species and a target species?**

The unit of biological explanation is therefore a directed source-target dyad, not a universal geographic field and not membership in a taxonomic group.

## Fresh hypothesis

For source species s and held-out target species t, define T_st as the transfer skill obtained when the spatial field is learned from source s and evaluated against target t after the inherited within-target IBD nuisance removal.

Before T_st is opened, define a response-blind relational similarity R_st from external ecological information.

Primary hypothesis:

    higher R_st -> higher T_st

conditional on predeclared geographic opportunity and broad lineage proximity.

A positive result licenses the statement that transferability is relational: ecological similarity predicts which source species carry reusable spatial information for a target species.

A null result licenses only that the frozen ecological relation did not predict transfer in the tested domain.

## Primary relational covariate: habitat / environmental-niche similarity

v0.1 uses **environmental-niche similarity** as the sole primary relation because it can be defined consistently across the full taxonomic panel.

The exact environmental variables, occurrence source, temporal filter, spatial thinning, background domain, minimum occurrence count, niche representation, overlap metric, missing-data rule, and transformations must be frozen before any genetic response for the fresh panel is opened.

Preferred implementation target:

- construct each species' environmental niche from external occurrence records;
- use one common predeclared environmental feature space for every species;
- define R_st with a symmetric continuous niche-overlap statistic;
- retain the continuous value; do not choose a high/low similarity threshold from genetic outcomes.

Candidate default for the next response-blind feasibility step is Schoener's D in a frozen environmental PCA space. This is **not yet authorized** until occurrence coverage and missingness are audited without genetic outcomes.

## Secondary relations

Host/resource similarity may be scientifically stronger for taxa where a defensible external interaction ontology exists, but it is not allowed to replace the primary relation after response opening.

It may enter only under a separately frozen secondary protocol with:
- taxon-specific applicability declared in advance;
- external provenance;
- a fixed similarity definition;
- a fixed missingness rule;
- no use of genetic response for coverage thresholds or ontology choices.

Taxonomic/phylogenetic similarity is a covariate/control, not the primary mechanism.

## Do not repeat the same-order source-pool contrast

The failed conditional v0.2 deployment changed source-pool composition and produced unacceptable private-null behavior on the confirmatory survivor geometry. Relational TTF therefore does not define the primary test as "high-similarity sources versus low-similarity sources."

Instead, all eligible source-target dyads are retained and the relation is tested against the dyad-level transfer matrix.

This preserves the scientific distinction between:
1. estimating T_st under a fixed transfer operator; and
2. explaining heterogeneity in T_st with pre-outcome R_st.

## Primary estimand

For eligible directed dyads (s,t):

    T_st = source-to-target post-IBD transfer skill
    R_st = frozen ecological similarity
    G_st = frozen geographic co-support
    P_st = frozen lineage/phylogenetic proximity

The conceptual primary model is

    T_st = a_t + g_s + beta_R R_st + beta_G G_st + beta_P P_st + error_st

where a_t captures target-specific predictability and g_s captures source-specific general usefulness.

Primary estimand: beta_R.

Primary directional hypothesis: beta_R > 0.

The final estimator need not be a parametric mixed model if synthetic qualification shows a more robust crossed-dyad statistic is required. The estimator and null generator must be frozen before empirical response opening.

## Dyadic dependence

Dyads are not independent because many observations share a source or a target. Ordinary edge/dyad bootstrap is forbidden.

Qualification must compare at least one predeclared crossed-source/target inference procedure that preserves source and target dependence. Candidate families are:
- crossed source-target random effects with calibrated inference;
- two-way source/target cluster bootstrap;
- source/target-label permutation that preserves the frozen dyadic geometry.

Selection among these is a response-blind development decision based only on synthetic Type-I/power and geometry, never empirical genetic performance.

## Geographic opportunity

Geographic co-support remains a nuisance/opportunity control, not the biological hypothesis.

R_st must not merely encode whether two species occur in the same places. The protocol must therefore freeze G_st and test ecological similarity beyond geographic opportunity. Environmental similarity and geographic co-occurrence must be reported separately.

## Fresh panels and firewall

Before panel construction, exclude every species whose empirical genetic response has contributed to:
- the completed 211-species unconditional result;
- any prior opened empirical genetic TTF analysis;
- any exploratory analysis that exposed pair-specific genetic transfer outcomes used to motivate this hypothesis.

Then create at least two species-disjoint sets under a deterministic response-blind hash ranking:

1. **development panel** — geometry, taxonomy, external ecology, and synthetic worlds only; empirical genetic responses remain closed;
2. **confirmatory panel** — untouched one-shot empirical test after all gates pass.

If an independent archive/source can supply a genuinely external confirmatory panel, prefer it over another slice of the same archive. If not, explicitly label same-archive confirmation as panel replication rather than source replication.

## Response-blind fresh-candidate census

The exact phylogatR source archive and prior-panel exclusions now yield a fresh 1,000-species geometry-eligible candidate set under the relational hash namespace. Reaching 1,000 eligible species required scanning 1,876 ranked candidates; 875 failed the inherited 12-locality geometry requirement and one lacked usable header-occurrence matches. The selected 1,000 have zero overlap with the prior 750 original-fresh and conditional-panel species.

The candidate pool spans 67 orders and 288 families. Lepidoptera is the largest order at 45%, followed by Diptera, Hymenoptera and Coleoptera. Because the intended primary claim is cross-taxon, taxonomic breadth is now a prospective ecological-feasibility guardrail rather than a post-result sensitivity choice:

- after external-ecology admissibility, no single order may exceed 50% of the fresh pool;
- at least four orders must each contribute at least 5%;
- failure closes the general cross-taxon route as `NOT_EVALUABLE_GENERAL_CROSS_TAXON_RELATIONAL_TTF`;
- a leave-largest-order-out analysis is secondary only and cannot rescue or overturn the primary decision.

The frozen candidate receipt is `benchmarks/frozen/relational_fresh_candidate_census_v0.1.json`. The full response-blind candidate table is retained outside the repository data source and can be regenerated with `scripts/freeze_relational_fresh_candidate_census.py`.

## Required pre-response gates

1. **Ecological-data feasibility** — enough species and dyads have externally defined R_st under the frozen missingness rule.
2. **Relation non-degeneracy** — R_st has sufficient response-blind variation within the eligible geography/lineage support.
3. **Geometry integrity** — no source/target dominance makes the estimand effectively a one-clade or one-region test.
4. **Private-null Type-I** — Wilson 95% upper bound <= 0.10 in every predeclared null cell.
5. **Relational positive-control power** — Wilson 95% lower bound >= 0.80 for the frozen minimally meaningful beta_R world.
6. **Fresh survivor requalification** — rerun exact gates after any character-mask filtering and before nucleotide identity / pairwise genetic distance opening.

Failure of a gate means NOT_EVALUABLE. It never authorizes tuning against empirical response.

## Interpretation branches

**Relational positive:** ecological source-target similarity predicts transfer skill beyond geography and broad lineage proximity. Spatial genetic information is conditionally reusable rather than universally transferable.

**Relational null with adequate power:** the frozen ecological relation does not explain source-target transfer heterogeneity in the tested panel.

**Gate failure:** no biological conclusion.

## Immediate next step

Run a response-blind feasibility census for environmental-niche similarity on a fresh candidate species universe. The census may inspect species identity, taxonomy, coordinates, occurrence availability, environmental coverage, and candidate R_st distributions. It must not compute or inspect fresh-panel genetic distances or T_st.

Only after that census should the exact niche metric and fresh development/confirmatory panel sizes be frozen.
