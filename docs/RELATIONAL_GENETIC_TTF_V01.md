# Relational genetic TTF v0.1 — fresh hypothesis contract

Status: **Study A empirical complete; host-resource relation null under qualified power. Study B environmental generalization remains response-blind.**

## Post-opening resolution — Study A

The prospectively frozen Lepidoptera host-resource study is complete. A response-independent mask-parser repair removed exactly three species that violated the pre-existing generic Phase-2 duplicate-aligned-header inadmissibility rule. The repair was frozen after the first identity-opening attempt stopped but before any source-target transfer score, `beta_R`, or relational decision existed. The repaired survivor geometry retained 420 mask-admissible species and an empirical union of 399 species, giving 203 source species, 196 target species, and 11,696 directed dyads.

The unchanged synthetic requalification passed on that repaired geometry: all four private-null cells had Wilson 95% rejection-rate upper bounds below 0.10, and the frozen small positive-control world had power 0.921 with Wilson 95% lower bound 0.903. The one authorized empirical continuation then resolved to:

- `beta_R = -0.00147409`;
- two-way source/target cluster-robust SE `= 0.00418434`;
- one-sided `p = 0.637688`;
- decision: `RELATIONAL_HOST_RESOURCE_NULL_WITH_QUALIFIED_POWER`.

Thus, within the frozen fresh Lepidoptera design, similarity in larval host-resource landscape footprints did **not** predict source-target transferability of post-IBD spatial genetic structure. This is not a rescue or reinterpretation of the earlier unconditional TTF result.

An exact reproduction audit confirmed all 11,696 repaired `T_st` values bit-for-bit against the authorized source vector and reproduced the same coefficient, SE and p-value. A later accelerated p-distance shortcut that altered floating-point summation/tie behavior is explicitly discarded for inference. A provenance audit also records that one estimator SHA string in the pre-response repaired-requalification receipt is unresolvable; rerunning the frozen synthetic worlds with the estimator blob explicitly bound by the v0.2 authorization reproduces every cell summary exactly.

Current authoritative receipts:
- `benchmarks/frozen/relational_host_resource_empirical_repaired_receipt_v0.2.json`;
- `benchmarks/frozen/relational_host_resource_empirical_reproduction_v0.2.json`;
- `benchmarks/frozen/relational_host_resource_provenance_audit_v0.2.json`.

No further host/resource retuning, alternative host metric selection, subgroup search, or result-selection rerun is authorized. Study B below remains a separately predeclared response-blind environmental-generalization study.


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

## Primary relational covariate resolved response-blind: host-resource similarity

The primary fresh test is now **Lepidoptera host-resource similarity**, not a global environmental PCA.

This choice was made before opening any fresh genetic response. A pre-existing host-resource design had already been constructed under an explicit zero-response firewall: sequence identity, pairwise genetic distance and transfer statistics were all unopened. It provides a direct biological relation rather than a generic similarity proxy.

The base host design contains 500 Lepidoptera. To make panel separation maximally conservative, 13 species that also appeared in the earlier conditional response-blind panel are dropped even though their nucleotide identity and genetic transfer outcomes were never opened. The resulting primary panel has:

- 487 species: 246 inherited source/train species and 241 inherited target/evaluation species;
- 226 evaluation targets with at least five geographically supported sources;
- 15,094 eligible directed source-target pairs;
- 47 families.

For each supported pair, the primary relation is

    R_host_st = Jaccard(primary-native larval host-resource WGSrpd3 footprint_s,
                        primary-native larval host-resource WGSrpd3 footprint_t)

It remains continuous. No high/low overlap threshold may be selected.

The response-blind distribution is strongly non-degenerate: the supported-pair Jaccard spans 0 to 1 with median approximately 0.350. Importantly, it is not simply a repackaging of geographic support: its response-blind Spearman association is about 0.043 with geographic coverage, -0.087 with source-target centroid distance, and near zero with edge/locality sampling-size ratios.

The exact primary contract is frozen in `docs/supporting/relational_host_resource_protocol_v0.1.json`, with its response-blind audit in `benchmarks/frozen/relational_host_resource_response_blind_audit_v0.1.json`.

## Habitat / environmental niche becomes external generalization

The 1,000-species cross-taxon candidate census remains useful, but it no longer gates the primary result. The panel mixes terrestrial, freshwater and marine animals, making a single terrestrial bioclimate PCA a weaker biological definition of habitat similarity.

The GBIF feasibility workflow therefore continues under a strict response firewall as a **future external-generalization study**. If a coherent cross-realm environmental representation can be frozen prospectively, it can ask whether relational transfer extends beyond Lepidoptera host resources. It cannot replace, rescue or reinterpret the host-resource primary test.

Schoener's D in an environmental PCA remains only a candidate for that successor and is not authorized by the host-resource study.

## Secondary relations

Within the primary Lepidoptera study, host breadth and host identity are secondary descriptors. Taxonomic family membership is a predeclared control. Any environmental relation remains a separately frozen external-generalization analysis and cannot be selected after genetic response opening.

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

Study A is closed. Continue only the separately predeclared Study B response-blind feasibility path: run the fresh 1,000-species GBIF occurrence/environment census without opening any genetic response for that panel. The census may inspect species identity, taxonomy, coordinates, occurrence availability, environmental coverage, and candidate ecological-relation distributions. It must not compute or inspect fresh-panel genetic distances or T_st.

If the frozen cross-taxon environmental representation is not ecologically coherent or fails its prospective feasibility/breadth gates, Study B is `NOT_EVALUABLE`; it must not be retuned into a favorable subset after the Study-A result.
