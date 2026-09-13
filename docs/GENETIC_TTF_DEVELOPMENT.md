# Genetic TTF development contract

Status: **development only; genetic outcomes not authorized for interpretation**.

This document defines the response-blind path for testing whether the geography of within-species genetic differentiation transfers to unseen evolutionary lineages. It does not declare the genetic interface qualified.

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

The main graph remains the core TTF species-local kNN graph. Development default: `k = 4`.

`src/ttf/genetic_geometry.py` freezes, without genetic outcomes, the original record count, exact unique localities, records per locality, kNN edge set, and for every edge the number of other edges sharing neither endpoint.

A coarse development pre-screen of **at least 12 exact unique localities** may be used to avoid obviously uninformative panels. It is not an empirical qualification criterion. Final EVALUABLE / NOT EVALUABLE status must be determined from the actual frozen kNN geometry by response-blind semi-synthetic calibration.

## 4. Biological IBD residualization

Pairwise genetic data have an additional dependence problem: edges that share a locality share biological information. Random edge cross-validation is therefore insufficient.

For each edge `e = (i, j)`, the genetic IBD nuisance fit excludes every other edge incident to `i` or `j`. On the remaining endpoint-disjoint edges it represents geographic and genetic distance only by training-set order, fits a rank-linear monotone IBD expectation, evaluates expected and observed order fractions for the held-out edge, and defines the residual as observed minus expected order fraction.

After every edge receives an out-of-pair residual, residuals are ranked within species. This is implemented in `src/ttf/genetic_ibd.py`.

The construction is invariant to strictly increasing rescalings of genetic or geographic distance and removes a noiseless strictly monotone IBD relation exactly.

This is a **biological estimand layer**, not the historical v0.3 geometry repair. The failed v0.3 fresh-external result remains failed and must not be reused as a qualification claim.

## 5. Relation to the current TTF inference lineage

Historical v0.3 controlled training-side edge-length rank but failed fresh external Type-I qualification on strong private spatial structure. v0.4 controlled Type I but lost too much power. The later profiled-private lineage reached fresh external Type-I validity, but deployment adequacy remained geometry specific: the 250-species RGFCA reserve passed Type I and failed the prespecified 80% power floor.

Therefore the genetic interface must not inherit qualification from any earlier panel. Its sequence is:

1. freeze genetic sampling geometry;
2. construct a genetic-specific semi-synthetic world including monotone IBD plus private and shared residual spatial transitions;
3. validate the cross-fit IBD response construction;
4. run the current profiled-private TTF validity machinery on the frozen genetic geometry;
5. require both Type-I control and prespecified positive-control power before any empirical genetic sharedness result is interpreted.

Failure at step 5 means **NOT EVALUABLE under this design**, not absence of shared phylogeographic structure.

## 6. Required synthetic arms before outcome opening

The genetic Gate-D surface must include at least:

- **IBD only:** strong monotone genetic distance with geographic distance, no place-specific residual field;
- **private residual:** every species has strong post-IBD spatial structure, but its location is species specific;
- **shared residual:** a prospectively defined fraction of species shares a common post-IBD transition field;
- **mixed IBD strength:** species differ in monotone IBD strength;
- **locality duplication:** multiple sequence records at the same locality do not inflate spatial information;
- **support mismatch:** inadequate train/evaluation geographic support yields NOT EVALUABLE rather than a biological negative.

Qualification thresholds must be frozen before observed genetic distances are opened.

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
