# TTF — Transferable Turnover Fields

TTF is a standalone, trait-agnostic methods repository for asking whether **within-species trait-transition structure learned from some species predicts transition structure in entirely unseen species**.

`fcp` remains an empirical flower-colour ecology project. TTF is the method-development, qualification, abstention, and transfer layer.

## Core estimand

TTF does **not** define cross-species sharedness as hotspot concentration in one pooled map.

> **Sharedness = out-of-species transferability.**

For species \(s\), TTF constructs a graph only among records of that species, measures edge-wise trait dissimilarity, rank-standardizes turnover within species, learns a field from training species, and scores that field on species held entirely outside training:

\[
C_s = \operatorname{Spearman}(\hat b_{se}, u_{se}), \qquad
T = |S_{\rm eval}|^{-1}\sum_s C_s.
\]

Earlier failed estimators and gates remain frozen. Later improvements never retroactively convert a failure into a pass.

## Core TTF qualification ledger

### v0.1 — trait-permutation null: immutable FAIL

The first prospective calibration showed that within-species trait permutation destroys strong private spatial structure and can therefore create anti-conservative inference. It remains a diagnostic exchangeability test, not the primary sharedness null.

### v0.2 — held-out-species inference: idealized PASS

The centered held-out-species bootstrap preserved private within-species spatial structure and passed the frozen high-precision idealized synthetic gate. This established the primary inferential object, but did not license arbitrary empirical geometries.

### v0.3–v0.10 — geometry and detectability development

These stages separated nuisance control from deployment detectability. A key geometry-only result was that increasing records from 20 to 100 with fixed `k=3` made the graph much more local: normalized edge scale fell to about **0.319** and shared-transition observability to about **0.335**. Keeping `k/n ≈ 0.15` retained both near **0.867**.

### v0.11 — fresh density-scaled confirmatory qualification: PASS

A prospectively frozen plant panel used **250 species × 100 records = 25,000 records**, a single density-scaled graph choice `k=15`, a 125/125 training/evaluation split, 500 observed worlds and 1,999 private-reference worlds.

Frozen result: `results/v11_fresh_density_scaled_qualification_v0.1.json`.

- private rejection at amplitudes 0.5 / 1 / 2 / 3: **0.046 / 0.070 / 0.020 / 0.050**;
- maximum Wilson 95% upper bound: **0.09580 ≤ 0.10**;
- shared amplitude-2 power: **0.998**;
- Wilson 95% lower bound: **0.98876 ≥ 0.80**.

The dense operator would require roughly 103.6 GiB, so the repository uses a chunked exact scorer checked against dense scoring to `1e-12`. This changes execution, not the estimand.

### v0.12 — actual observation-support gate: NOT EVALUABLE

The same frozen 25,000 photos were next passed through an independently frozen, location-blind measurement pipeline **before** TTF was allowed to read colour vectors, pairwise colour distances, edge turnover, or an empirical TTF statistic.

The prospective requirement was that all 250 species retain at least 40 evaluable photographs. Terminal result:

- species meeting the minimum: **189 / 250**;
- median evaluable photographs/species: **51**;
- minimum: **15**;
- retained actual geometry: **0 species / 0 records**;
- `colour_vector_values_read_by_ttf = false`;
- `pairwise_colour_distances_computed = false`;
- `ttf_empirical_statistic_computed = false`.

Therefore v0.12 is a **measurement-support / observational-admissibility failure**, not an ecological null and not evidence that shared flower-colour transitions are absent. The correct endpoint is **NOT EVALUABLE / ABSTAIN**.

Terminal ledger: `benchmarks/frozen/v12_actual_geometry_source.json`.

## Paired-state extension: TTF-M and TTF-C

TTF also distinguishes two paired-state transfer estimands:

- **TTF-M:** transfer turnover of a predeclared mismatch magnitude `M_s(x)=m(A_s(x),B_s(x))`;
- **TTF-C:** transfer turnover of a predeclared relation state `R_s(x)=r(A_s(x),B_s(x))`.

The direction-free idealized qualification passed for both. Under shared relation rotation with constant mismatch magnitude, TTF-M rejected `0/500` while TTF-C rejected `500/500`, demonstrating that magnitude-turnover and relation-turnover are distinct targets.

### Q4 directionality

Q4 v0.1 is an immutable FAIL because its supposedly private worlds shared a central front-density zone. Q4.1 replaced that world definition prospectively and PASSed: shared BREAK and RECOUPLE were detected `500/500`; neutral rotation and sign-conflict received no directional label; strict-private relation/breakdown rejection was `24/500 = 0.048`, Wilson upper `0.07043`.

`breakdown` is only a transferable increase in a predeclared mismatch coordinate across a predeclared orientation. It is not, by itself, causality, fitness loss, or mechanism failure.

### Heterogeneous geometry: terminal current ceiling

The heterogeneous-geometry successor sequence is now closed:

```text
Q5 FAIL
  -> non-qualifying factor diagnosis
  -> Q5.1 training-only edge-length control FAIL
  -> Q5.2 response-blind geometry transport: candidate_selected = null
```

Q5/Q5.1 establish a stable asymmetry on the declared finite synthetic families:

- **TTF-M: PASS**;
- **TTF-C: FAIL**;
- **joint extension: FAIL**.

Q5.1 improved matched shared-mismatch power but still failed private type-I control and shifted shared-mismatch power. Q5.2 then prospectively tested three field-mass transports — `pooled_inverse_density`, `self_inverse_density`, and `target_density_ratio` — on fresh worlds. None met all development margins, so no formal Q5.2 run was permitted.

Frozen Q5.2 result: `results/ttfm_q5_2_geometry_transport_development_v0.1.json`.

This means **heterogeneous-geometry TTF-C remains unqualified**. The failed family may not be rescued by post-outcome threshold, graph, bandwidth, clipping, response-residualization, null, or candidate-combination tuning. Any future successor requires a newly named mechanism diagnosis and a fresh development family.

## Why the evidence layers stay separate

TTF explicitly separates:

```text
method validity on declared worlds
        !=
deployment-geometry detectability
        !=
actual observational admissibility
        !=
empirical biological conclusion
```

A method can pass a synthetic qualification and still be unusable for a specific empirical dataset because measurement support failed. Likewise, one transfer estimand can survive heterogeneous geometry while another does not.

## PAYOFF mechanistic spatial projection

TTF can also serve as the **cross-species spatial-projection layer** for an independently specified architecture-payoff mechanism.

```text
Campanula microdonta
individual-system anchor
        |
        v
PAYOFF local architecture margin
        |
        v
spatial payoff projection
        |
        v
TTF held-out cross-species transfer
```

The generic predictor-space API accepts arbitrary named feature spaces, so a prospectively specified PAYOFF-derived predictor `P` can be compared with geography/environment baselines on the same held-out species split. For example, `P|GE = T_GEP - T_GE` asks whether a mechanism-derived projection adds transfer skill beyond the declared geography/environment space.

Held-out trait outcomes may not be used to tune the PAYOFF feature. The Izu *Campanula microdonta* populations are one deeply resolved anchor species, not multiple independent TTF species.

See `docs/PAYOFF_SPATIAL_PROJECTION_HANDOFF.md` for the feature contract, anti-leakage rules, paired increments, directional caveat, and claim boundary.

## Implemented

- within-species-only kNN graphs;
- generic scalar/vector trait dissimilarity;
- within-species rank standardization;
- species-equal opportunity-corrected kernel fields;
- exact chunked scoring for large frozen geometries;
- species-disjoint transfer statistics and held-out-system bootstrap inference;
- frozen synthetic and empirical-geometry qualification workflows;
- TTF-M / TTF-C paired-state estimands and bounded directional labels;
- generic held-out predictor-space competition with paired increments;
- explicit observation-support abstention gates.

See `docs/METHOD_SPEC.md`, `docs/QUALIFICATION_PROTOCOL.md`, and `docs/TTFM_METHOD_SPEC.md`.

## Install and test

```bash
python -m pip install -e ".[test]"
pytest -q
```

## Current claim ceiling

TTF currently licenses bounded method statements only.

- Core TTF v0.11 is qualified on its declared synthetic/frozen-geometry design.
- The actual flower-colour route stopped at v0.12 as **NOT EVALUABLE before any empirical TTF statistic was opened**.
- TTF-M and TTF-C are both qualified on their balanced idealized family, and Q4.1 passes its strict-private directional family.
- Under the declared heterogeneous Q5/Q5.1 families, **TTF-M passes while TTF-C remains unqualified**; Q5.2 selected no successor candidate.
- PAYOFF projection is an external mechanistic predictor contract, not evidence that the mechanism has already been empirically validated.

No empirical flower-colour sharedness/no-sharedness, pollinator, island, fitness, causal mismatch, or universal robustness conclusion is licensed by this repository state.