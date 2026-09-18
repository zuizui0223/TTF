# Genetic TTF development contract

Status: **historical response-blind development contract — preserved after empirical completion**.

**Post-opening resolution (2026-09-18):** the independent fresh phylogatR route completed under this development lineage. The exact survivor geometry passed cross-species Gate-D and self-detectability qualification, and the single authorized Phase-4 opening resolved to `LINEAGE_CONDITIONED_SPATIAL_STRUCTURE_WITHIN_TESTED_DOMAIN` (`place_beyond_ibd T=0.0325655`, profiled-private `p=0.739261`; self `S=-0.226098`, upper-tail `p=0.00899101`). See `benchmarks/frozen/genetic_phylogatr_phase4_empirical_handoff_v0.1.json` and `manuscript/genetic_ttf_flagship_v0.2.md`. The development-only rules below are retained as the audit trail and must not be read as the current empirical state.

This document defines the response-blind development path for testing whether the geography of within-species genetic differentiation transfers to unseen evolutionary lineages.

## 1. Scientific target

The primary question is:

> Does geographic location predict within-species genetic differentiation in an unseen species beyond ordinary isolation by distance (IBD)?

Two estimands are separated:

1. **total genetic transfer** — transferability of within-species genetic-distance rank before removing biological IBD;
2. **place beyond IBD** — transferability of genetic differentiation remaining after a prospectively fixed, cross-fit within-species IBD expectation is removed.

A positive total result with no residual result is evidence for generic distance / geometry, not for a transferable place-specific field.

## 2. Development data firewall

The current public development source is the sampling geometry in `skdecker/PhylogeographicBreaks`, frozen for this development pass at commit `11534ba70fbe3705e667651030edef5b849fa949`.

Allowed before the geometry contract is frozen: species/locus/accession identifiers, coordinates, exact coordinate duplication, record counts, kNN graph geometry and edge lengths, and response-blind support overlap.

Forbidden before authorization to open genetic outcomes: pairwise sequence distances, alignment-derived divergence summaries, Monmonier break results, break presence/absence, or tuning eligibility thresholds using any genetic outcome.

GitHub code-search census found 230 Aves and 165 Chiroptera candidate `*_occurfiltered.csv` paths. These are candidate occurrence files, not yet a frozen analysis manifest.

The Decker data are a **development/pilot source**, not the preferred confirmatory source, because their upstream sequence filtering was designed for a different phylogeographic analysis. A confirmatory analysis should use a fresh phylogatR extraction after all genetic-distance and qualification rules are frozen.

## 3. Response-blind locality geometry

Sequence count is not an information criterion for TTF. Multiple sequences may share one coordinate, so records are first collapsed by **exact coordinate pair**. Near-but-nonidentical points are not merged unless a separate radius rule is frozen before genetic outcomes are opened.

The genetic main graph now inherits the **qualified core TTF v0.11 density-scaled architecture**, rather than the historical fixed `k=4` development default. After exact locality collapse,

`k_s = min(n_s - 1, max(2, floor(0.15 * n_s + 0.5)))`.

This rule is determined by locality geometry only. It follows the v0.10 mechanism diagnosis and v0.11 fresh confirmatory PASS, where fixed small `k` became too local as sampling density increased while `k/n ~= 0.15` preserved transition observability.

`src/ttf/genetic_geometry.py` freezes, without genetic outcomes, the original record count, exact unique localities, records per locality, selected density-scaled `k`, kNN edge set, and for every edge the number of other edges sharing neither endpoint.

A coarse development pre-screen of **at least 12 exact unique localities** may be used to avoid obviously uninformative panels. It is not an empirical qualification criterion. Final EVALUABLE / NOT EVALUABLE status must be determined from the actual frozen density-scaled geometry by response-blind semi-synthetic calibration.

## 4. Biological IBD residualization

Pairwise genetic data have an additional dependence problem: edges that share a locality share biological information. Random edge cross-validation is therefore insufficient.

For each edge `e = (i, j)`, the genetic IBD nuisance fit excludes every other edge incident to `i` or `j`. On the remaining endpoint-disjoint edges it represents geographic and genetic distance only by training-set order, fits a rank-linear monotone IBD expectation, evaluates expected and observed order fractions for the held-out edge, and defines the residual as observed minus expected order fraction.

After every edge receives an out-of-pair residual, residuals are ranked within species. This is implemented in `src/ttf/genetic_ibd.py`.

The construction is invariant to strictly increasing rescalings of genetic or geographic distance and removes a noiseless strictly monotone IBD relation exactly.

This is a **biological estimand layer**, not a geometry repair. Core TTF geometry safeguards are applied only after the biological IBD expectation has been removed.

## 5. Relation to the current TTF inference lineage

The relevant core lineage is now v0.11, not v0.3 or v0.8.

- v0.3/v0.4 remain immutable development failures.
- v0.8 established that profiled-private inference can control private-structure Type I on fresh external panels, but deployment power remained geometry specific.
- v0.10 identified the edge-scale mechanism: fixed small `k` becomes increasingly local as record density rises.
- **v0.11 passed** a fresh 250-species x 100-record confirmatory qualification using the single frozen density-scaled choice `k=15`, i.e. `k/n=0.15`, together with training-only edge-length orthogonalization and profiled-private inference.
- v0.12 stopped the flower-colour empirical route at an observation-support gate; that measurement failure is not a failure of the v0.11 core method.

The genetic interface therefore inherits the *architecture* of v0.11, but not its qualification. Its sequence is:

1. freeze exact-locality sampling geometry and density-scaled graphs;
2. require sufficient endpoint-disjoint edges for leave-two-localities-out IBD fitting;
3. generate genetic-specific semi-synthetic worlds with monotone IBD plus private/shared residual differentiation;
4. construct post-IBD residual turnover;
5. on training species only, apply the inherited edge-length-rank orthogonalization;
6. use the inherited training-only private-strength profile and profiled-private null;
7. require dataset-specific Type-I control **and** positive-control power before any observed genetic transfer statistic is interpreted.

`src/ttf/genetic_simulate.py` implements the synthetic genetic worlds from latent locality states, so pairwise distances share endpoints rather than being independent edge draws. `src/ttf/genetic_gate.py` implements the combined post-IBD v0.11 scoring and private-reference machinery.

Failure at step 7 means **NOT EVALUABLE under this design**, not absence of shared phylogeographic structure.

## 6. Required synthetic arms before outcome opening

The genetic Gate-D surface must include at least:

- **IBD only:** strong monotone genetic distance with geographic distance, no place-specific residual field;
- **private residual:** every species has strong post-IBD spatial structure, but its location is species specific;
- **shared residual:** a prospectively defined fraction of species shares a common post-IBD transition field;
- **mixed IBD strength:** species differ in monotone IBD strength;
- **locality duplication:** multiple sequence records at the same locality do not inflate spatial information;
- **support mismatch:** inadequate train/evaluation geographic support yields NOT EVALUABLE rather than a biological negative.

The initial synthetic generator represents each locality by a latent Euclidean genetic state. Its geographic dimensions produce exact monotone IBD when residual amplitude and noise are zero; a separate transition dimension creates shared/private residual differentiation; locality-level noise preserves pairwise endpoint dependence. Thus the IBD-only arm can test the whole response-construction pipeline rather than merely unit-testing a regression helper.

Qualification thresholds must be frozen before observed genetic distances are opened. The intended high-precision convention is the same as core TTF: private-null Wilson 95% upper bounds must remain below the frozen Type-I ceiling, and the shared moderate-signal Wilson 95% lower bound must exceed the frozen power floor.

## 7. Empirical hierarchy after qualification

If and only if the genetic geometry qualifies:

1. random species-disjoint train/evaluation transfer;
2. phylogenetic-block transfer;
3. Aves -> Chiroptera;
4. Chiroptera -> Aves.

Cross-class tests are the strongest guard against explaining a positive result only by close relatives sharing demographic history. They remain conditional on spatial support overlap and adequate power in each direction.

## 8. Claim boundary

Until genetic-specific qualification passes, this branch supports only:

> We have a response-blind, endpoint-safe interface for asking whether residual phylogeographic differentiation is spatially transferable across species.

It does **not** yet support an empirical claim that such transferability exists.
