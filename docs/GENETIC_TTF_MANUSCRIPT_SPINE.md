# Genetic TTF flagship manuscript spine

Status: **outcome-blind pre-data manuscript contract**.

This document fixes the biological story, estimands, section order, figure logic, and result-interpretation branches for the genetic TTF flagship **before** a fresh confirmatory nucleotide identity or empirical TTF statistic is opened.

It is not an empirical result and it does not authorize any new data opening. The canonical empirical opening state remains `benchmarks/frozen/genetic_empirical_opening_state_v0.3.json`.

---

## 1. Working paper identity

### Working title

**A transferable geography of intraspecific genetic differentiation**

Alternative subtitle-style title if a methods-oriented framing is needed later:

**Testing whether phylogeographic structure transfers to unseen species beyond isolation by distance**

### One-sentence biological question

> Does geographic location predict where within-species genetic differentiation occurs in evolutionary lineages that were not used to learn that geography, beyond ordinary isolation by distance?

### Conceptual contrast

The paper is not asking whether several co-distributed species show similar historical breaks after the fact.

It asks a stronger predictive question:

> Can a spatial turnover field learned from some species predict differentiation in entirely unseen species?

The distinction is essential. Sharedness is defined operationally as **out-of-species predictive transfer**, not as pooled hotspot concentration, repeated visual overlap, or retrospective concordance of named breaks.

---

## 2. Biological motivation

Population differentiation is spatially structured by dispersal limitation, barriers, climatic history, topography, habitat transitions, and lineage-specific demographic history. Comparative phylogeography has repeatedly shown that co-distributed species may exhibit concordant spatial structure, but concordance alone does not identify how much information belongs to geography itself rather than to the histories of the particular lineages being compared.

A location-based component should satisfy a stronger criterion: information learned from one set of species should predict differentiation in species that were completely absent from training.

Genetic distance creates an additional nuisance. Ordinary isolation by distance can generate broad spatial predictability even when no place-specific transition field exists. Therefore the flagship estimand is not generic genetic-distance transfer, but **place beyond IBD**: transferable spatial differentiation after an endpoint-safe within-species IBD expectation has been removed.

This converts a familiar comparative-phylogeographic idea into a prospective predictive test with an explicit abstention state when the sampling geometry cannot support the inference.

---

## 3. Fixed estimands

Two estimands remain separate throughout the manuscript.

### 3.1 Total genetic transfer

Transferability of within-species genetic-distance rank before removing biological isolation by distance.

This estimand can detect generic distance/geometry effects and is therefore secondary for the main biological claim.

### 3.2 Place beyond IBD — primary estimand

Transferability of genetic differentiation remaining after a prospectively fixed, leave-two-localities-out within-species IBD expectation is removed.

Interpretation is frozen as follows:

- total positive + residual positive: evidence for a transferable place component beyond generic IBD;
- total positive + residual null: evidence compatible with generic distance/geometry but not a transferable place-specific field;
- residual non-significant + qualified within-species self-detectability positive: evidence consistent with lineage-conditioned spatial structure within the tested domain;
- residual non-significant without qualified self support: `NOT_EVALUABLE_FOR_LINEAGE_CONDITIONING`;
- any failed geometry, admissibility, Type-I, or power gate: `NOT_EVALUABLE`, never a biological null.

The manuscript must not collapse these states into a single positive/negative dichotomy.

---

## 4. What is new relative to ordinary comparative phylogeography

The manuscript's novelty claim should be built around four differences.

1. **Prediction rather than retrospective concordance.** Species used for evaluation are entirely excluded from field learning.
2. **A spatial field rather than a species-level break label.** The response is edge-level within-species differentiation over geography, not simply whether a species has a detected phylogeographic break.
3. **Place beyond IBD.** Generic monotone distance dependence is removed using an endpoint-safe nuisance fit before testing transfer.
4. **Qualification before interpretation.** The exact empirical sampling geometry must demonstrate private-null Type-I control and positive-control power before a biological result can be interpreted.

Decker et al. provide the nearest development comparison because they predict phylogeographic-break occurrence from species characteristics, but they do not test whether a spatial differentiation field learned from some species predicts the location of differentiation in unseen species. Their bird/bat data remain development-only here.

---

## 5. Fixed data hierarchy

### 5.1 Development source

The Decker bird/bat panel is retained only as a development geometry.

Its frozen cross-species Gate-D result is `NOT_EVALUABLE` because the private A3 Wilson 95% upper bound (`0.10033475332223055`) narrowly exceeds the frozen `0.10` Type-I ceiling. The shared-A2 power gate passed (`Wilson lower = 0.8355052092686175`).

The same development geometry passes the separate within-species self-detectability qualification (`null Wilson upper = 0.053771365501634506`; `private-A2 Wilson lower = 0.9923756595384479`).

Therefore the development failure is informative about deployment inference, not evidence that within-species spatial genetic structure is absent.

Decker nucleotide identity, pairwise genetic distances, Monmonier outcomes, and break labels remain closed.

### 5.2 Confirmatory source

The confirmatory dataset must be a genuinely fresh authenticated phylogatR archive independent of the Decker development files.

The frozen route is:

- Animalia;
- exclude Aves;
- exclude Chiroptera;
- exclude the exact frozen 221 Decker development species;
- COI/COX1 marker family under the frozen alias rule;
- one geometry-selected panel per species;
- at least 12 exact unique localities;
- density-scaled within-species graph with `k/n = 0.15`;
- at least 5 endpoint-disjoint IBD training edges;
- maximum 250 species;
- minimum 150 species;
- deterministic 50/50 species split;
- no post-outcome support dropping, marker switching, rewiring, or resplitting.

The current blocker is source acquisition only. Anonymous phylogatR access, BOLD contingency access, and current public-GitHub raw-source substitution have already been closed or placed on hold under frozen response-blind preflights.

---

## 6. Response construction

### 6.1 Locality is the node

The biological node is an exact unique geographic locality, not a sequence record. Multiple sequences at one coordinate do not create additional spatial nodes.

### 6.2 Genetic edge response

For each frozen locality pair, sequence-pair genetic distance is the uncorrected p-distance over jointly canonical A/C/G/T columns, with the frozen comparable-site threshold. The locality-pair response is the arithmetic mean over all valid cross-locality sequence pairs.

If a frozen edge has no valid sequence pair, the species becomes empirical-outcome `NOT_EVALUABLE`; the graph is not rewired.

### 6.3 Endpoint-safe IBD removal

For a target edge `(i,j)`, every nuisance-fit edge incident to locality `i` or `j` is excluded. The IBD expectation is learned only from endpoint-disjoint within-species edges.

This prevents the held-out edge from borrowing biological information through a shared locality.

After all edges receive cross-fit residuals, residual turnover is rank-standardized within species.

---

## 7. Transfer model and inference

The inherited core architecture is the qualified v0.11 TTF design, but the genetic dataset must qualify independently.

The fixed sequence is:

1. construct frozen exact-locality density-scaled graphs;
2. compute endpoint-safe post-IBD residual turnover;
3. learn the spatial turnover field from training species only;
4. apply the inherited training-only edge-length-rank control;
5. predict residual turnover on species excluded entirely from training;
6. summarize each evaluation species by Spearman correlation between predicted and observed residual turnover;
7. average equally across evaluation species;
8. evaluate the statistic against the frozen profiled-private reference distribution.

The primary statistic is therefore a species-disjoint transfer score, not a pooled edge correlation.

---

## 8. Prospective opening sequence

The manuscript should report the data-opening sequence explicitly because it is part of the inferential design.

### Phase 1 — geometry only

Allowed: provenance, accession/header information, coordinates, geometry, file hashes.

Forbidden: nucleotide identity.

Pass requirement: `FROZEN_RESPONSE_BLIND_PHASE1_GEOMETRY` with at least 150 species.

### Phase 2 — canonical-valid mask only

Allowed: A/C/G/T-versus-noncanonical validity masks.

Forbidden: nucleotide identity.

Every frozen edge must retain at least one valid cross-locality sequence pair under the frozen comparable-site rule. No graph repair is permitted.

### Phase 3 — exact-geometry synthetic qualification

The full dataset-specific Gate-D is run on the exact Phase-2 survivor geometry.

Frozen high-precision requirements:

- private-null Wilson 95% upper bound `<= 0.10`;
- shared-A2 Wilson 95% lower bound `>= 0.80`.

Failure is `NOT_EVALUABLE` and nucleotide identity remains closed.

A fresh within-species self-detectability qualification is also run on the exact fresh geometry. It is required only for interpreting a later qualified cross-species non-significant result as lineage-conditioned.

### Phase 4 — nucleotide identity and one empirical test

Nucleotide identity can be opened only after complete Phase-3 PASS on the exact survivor geometry.

The empirical response is computed in memory under the frozen p-distance contract and the primary profiled-private TTF is evaluated once.

---

## 9. Outcome-independent manuscript architecture

### Introduction

**Paragraph 1 — biological problem.** Geographic structure in within-species genetic variation reflects both repeatable properties of place and lineage-specific history.

**Paragraph 2 — limitation of current comparative logic.** Repeated breaks and concordance are informative but retrospective; they do not test whether geography learned from one lineage predicts another.

**Paragraph 3 — predictive criterion.** Define a transferable geographic component as out-of-species prediction, and separate generic IBD from place-specific transfer.

**Paragraph 4 — study design.** Introduce endpoint-safe genetic TTF, prospective synthetic qualification on the exact empirical geometry, and a fresh unseen-lineage confirmatory test.

### Methods

1. Predictive definition of shared phylogeographic structure.
2. Development-versus-confirmatory data separation.
3. Exact-locality graph construction.
4. Genetic distance and admissibility contract.
5. Endpoint-safe IBD residualization.
6. Cross-species turnover field.
7. Profiled-private inference.
8. Exact-geometry synthetic qualification.
9. Within-species self-detectability.
10. Phased outcome firewall and stopping rules.

### Results

The first results section is fixed regardless of the eventual empirical outcome:

**R1. A generic TTF architecture can be adapted to pairwise genetic data only after endpoint-safe IBD control.**

**R2. Development geometry demonstrates why dataset-specific qualification is necessary.** Report Decker Gate-D `NOT_EVALUABLE` and self-detectability `PASS` without opening its empirical genetic outcomes.

**R3. Fresh confirmatory geometry and admissibility.** Report Phase-1/2 panel size, locality/edge distributions, support, and any pre-outcome exclusions.

**R4. Exact-geometry qualification.** Report all frozen private-null cells, shared-A2 power, and fresh self-detectability qualification.

**R5. Empirical out-of-species transfer.** This section is populated only if Phase 3 authorizes Phase 4.

No earlier section may be rewritten to make R5 appear inevitable after its outcome is known.

---

## 10. Frozen empirical result branches

### Branch A — transferable place component

Condition:

- fresh cross-species qualification PASS;
- primary `place_beyond_ibd` empirical `p <= 0.05`.

Licensed headline:

> Geographic location contains transferable information about within-species genetic differentiation that generalizes to unseen species beyond ordinary isolation by distance.

Not licensed without additional evidence:

- a specific causal barrier;
- common demographic history;
- identical historical mechanism among species;
- a universal global field.

### Branch B — lineage-conditioned spatial structure

Condition:

- fresh cross-species qualification PASS;
- primary empirical `p > 0.05`;
- fresh self-detectability qualification PASS;
- empirical within-species self statistic significant under its frozen test.

Licensed headline:

> Within-species differentiation is spatially predictable within lineages but does not transfer detectably across unseen species under the tested geographic field, consistent with lineage-conditioned spatial structure in this domain.

This is not a proof that geography never transfers in other clades, markers, scales, or regions.

### Branch C — cross-species non-significant but lineage conditioning not diagnosable

Condition:

- cross-species qualification PASS;
- primary empirical `p > 0.05`;
- self qualification or empirical self support absent.

Required endpoint:

`NOT_EVALUABLE_FOR_LINEAGE_CONDITIONING`.

The paper may report a qualified non-significant cross-species transfer test, but it must not convert that result into evidence for lineage-specific history.

### Branch D — confirmatory design not evaluable

Condition: Phase 1, Phase 2, Phase 3 Type-I, or Phase 3 power gate fails.

Required endpoint:

`NOT_EVALUABLE` / abstain.

This branch supports a methodological/deployment conclusion about the sampled geometry, never a biological claim of no transferable phylogeographic structure.

---

## 11. Figure plan

### Figure 1 — What would it mean for geography to transfer across lineages?

Conceptual map with training species and entirely unseen evaluation species. Contrast retrospective overlapping breaks with prospective prediction of held-out differentiation.

### Figure 2 — Endpoint-safe genetic TTF

Localities -> density-scaled edges -> pairwise genetic response -> leave-two-localities-out IBD expectation -> residual turnover -> training-species field -> unseen-species prediction.

The visual should make shared-endpoint leakage impossible to miss.

### Figure 3 — Qualification before biology

Development Decker geometry as the cautionary example: self-detectability can pass while cross-species Type-I qualification narrowly fails. Then show the fresh confirmatory opening sequence Phase 1 -> 2 -> 3 -> 4.

### Figure 4 — Fresh confirmatory field and transfer

Only after authorization. Show training field, held-out species coverage/support, per-species transfer coefficients, and the primary aggregate statistic against its frozen private reference.

### Figure 5 — Place versus lineage interpretation

Decision diagram using cross-species transfer and within-species self-detectability to distinguish transferable place component, lineage-conditioned structure, and abstention.

If the eventual empirical branch is simple enough, Figures 4 and 5 may be combined; the inferential states themselves must remain explicit.

---

## 12. Supplementary structure

- Table S1: frozen protocols, receipts, commits, and opening states.
- Table S2: Phase-1 species/panel geometry ledger.
- Table S3: Phase-2 admissibility ledger without nucleotide identities.
- Table S4: all Phase-3 private-null and shared-positive-control cells.
- Table S5: fresh self-detectability qualification.
- Table S6: per-species empirical transfer coefficients if Phase 4 is authorized.
- Figure S1: graph-density and endpoint-disjoint-edge distributions.
- Figure S2: support-overlap diagnostics.
- Figure S3: total genetic transfer versus place-beyond-IBD estimands.
- Figure S4: sensitivity outputs that were prospectively authorized; no post-outcome rescue analyses.

---

## 13. Abstract template fixed before outcome

**Background:** Comparative phylogeography asks whether species share spatial genetic structure, but retrospective concordance cannot by itself determine whether predictive information belongs to geography or to the histories of the lineages compared.

**Approach:** We define shared geographic structure as out-of-species transferability. We learn a spatial field of within-species genetic differentiation from training species and test whether it predicts differentiation in unseen species after endpoint-safe removal of isolation by distance. The exact empirical geometry must pass prospective Type-I and power qualification before nucleotide identities can be opened for the confirmatory test.

**Result sentence:** populate only from Branch A, B, C, or D above.

**Conclusion sentence:** must remain at the same inferential level as the selected branch and may not introduce a stronger causal mechanism than was tested.

---

## 14. Current claim ceiling before fresh archive acquisition

At the current repository state the manuscript may state only that:

1. a genetic TTF interface has been implemented for testing out-of-species transfer beyond IBD;
2. the construction uses endpoint-safe within-species IBD removal and a species-disjoint transfer estimand;
3. the Decker development geometry is `NOT_EVALUABLE` for cross-species inference under the frozen Gate-D despite passing within-species self-detectability;
4. a prospective fresh phylogatR confirmatory route is frozen through Phase 4;
5. public-source alternatives have been closed or placed on hold under response-blind preflights;
6. no empirical genetic transferability conclusion currently exists.

The manuscript must not imply that a transferable geography has already been observed.

---

## 15. Stop rules for manuscript development

Before the fresh empirical outcome is opened:

- do not write a result-specific title;
- do not select a biological mechanism that only becomes attractive after seeing the transfer map;
- do not redefine the primary estimand;
- do not change graph density, bandwidth, IBD rule, private-null construction, alpha, or qualification thresholds to improve the eventual observed result;
- do not exclude inconvenient held-out species after outcome opening unless a prospectively frozen admissibility rule already requires it;
- do not convert `NOT_EVALUABLE` into a negative biological result;
- do not use the development Decker outcomes as a substitute empirical result.

The next legitimate empirical change to this manuscript begins only after an untouched fresh phylogatR archive passes the frozen response-blind Phase-1 intake.
