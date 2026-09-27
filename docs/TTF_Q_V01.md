# TTF-Q v0.1 — response-blind qualification of dyadic ecological predictors

Status: **method-development branch, separate from the permanently closed relational v0.3 program**.

## Why this exists

Relational TTF v0.3 ended in `CLOSED_NO_EVALUABLE_TEST`. Study B and Study C were not biological nulls: both stopped on response-blind external-data/transport requirements before nucleotide identity, pairwise genetic distances, `T_st`, or the focal genetic coefficients were opened.

TTF-Q does **not** reopen those slots. It extracts the reusable methodological problem that the frozen B/C pipeline exposed:

> Before a source-target ecological relation is tested against an outcome, can we determine whether that relation is measured stably, contains information distinct from known structure, is distributed over enough source/target support, and could detect biologically meaningful effects under the intended dependence model?

The output is an **evaluable envelope**, not `PASS`, `FAIL`, `significant`, or `null`.

## Core object

For directed source-target dyads `(s,t)`, let `R_st` be a focal ecological relation. Examples are:

- B: Schoener's D for current environmental niches, `R_env`.
- C: Schoener's D for historical climate-displacement distributions, `R_hist`.

The intended model may also contain geographic opportunity, lineage indicators, locality-count imbalance, and exact source and target fixed effects.

TTF-Q characterizes four response-blind properties.

## 1. Relation repeatability

When the relation can be rebuilt under occurrence/bootstrap resampling, TTF-Q receives a matrix whose rows are resampling replicates and columns are the same dyads.

It reports the balanced one-way random-effects repeatability ICC plus within-dyad uncertainty. This asks whether a dyad ranked as highly similar remains highly similar when the external observations used to estimate the relation are perturbed.

This is a measurement diagnostic. It does not use or predict the biological response.

## 2. Non-redundant relational information

The focal relation is first put on the exact scale used by the intended model. Source and target fixed effects are then absorbed exactly. The remaining focal relation is projected on the intended nuisance controls.

For focal vector `r`, two-way-FE residual `r_FE`, and control-residualized vector `r_unique`, TTF-Q reports:

- fraction retained after source/target identity: `||r_FE||^2 / ||r||^2`;
- fraction unique after controls: `||r_unique||^2 / ||r_FE||^2`;
- total unique fraction: `||r_unique||^2 / ||r||^2`;
- a VIF-like quantity `1 / unique_fraction_after_FE`.

### Exact information-loss decomposition

The three fractions are not unrelated diagnostics. By construction they satisfy the exact identity

```
total_unique_fraction
= source_target_FE_retained_fraction
  × control_unique_fraction_after_FE
```

or, equivalently,

```
||r_unique||² / ||r||²
= (||r_FE||² / ||r||²)
  (||r_unique||² / ||r_FE||²).
```

This separates two different reasons a relational predictor can become uninformative:

1. **endpoint-identity loss** — most apparent relation variation is actually source or target identity;
2. **nuisance-redundancy loss** — the within-endpoint relation is largely reproduced by geography, lineage, current climate, or other declared controls.

The decomposition is exact for the declared linear projection and therefore turns “not enough independent relation information” into a localized diagnosis rather than a binary failure label.

For Study C this directly answers the scientifically useful question that survives the closed confirmatory program:

> How much historical-climate similarity exists that cannot be reduced to current-climate similarity, geography, lineage, locality imbalance, or source/target identity?

That is a property of the ecological relation itself and can be studied without genetic outcomes.

## 3. Dyadic signal support

A large dyad count can hide concentration in a handful of source or target species. TTF-Q therefore uses squared non-redundant focal-relation signal and reports inverse-Herfindahl effective counts across:

- dyads;
- source endpoints;
- target endpoints;

plus the maximum source and target signal shares.

These are **signal-concentration diagnostics**, not claims that dyads are independent and not replacements for cluster-robust inference.

## 4. Detectability surface

The frozen B/C qualification code reduced a rich synthetic experiment to a binary gate at one effect size. TTF-Q retains the same intended two-way source/target dependence-aware estimator but reports the full surface.

The v0.1 case-study grid is:

- standardized focal effect: `0, .01, .02, .03, .05, .08, .10`;
- private-heterogeneity amplitude: `A = 0, 1, 2, 3`;
- source and target intercept SD: `0.18`;
- dyad noise SD: `0.32`;
- private random-slope SD: `0.12 * A`;
- response transform: `tanh`;
- intended inference alpha: `0.025`.

For every cell TTF-Q reports rejection rate and its Wilson 95% interval. It then reports, for each `A`, the smallest **grid** effect whose Wilson lower bound reaches the target power of 0.80.

This converts the old question

> Did beta=.03 pass the gate?

into the more useful question

> Under this actual dyadic geometry, what effect sizes are detectable as source/target heterogeneity worsens?

No ecological effect is inferred from the synthetic worlds.

## B and C as case studies

### B — current environmental similarity

The focal relation is `z_R_env`. Controls are geographic coverage, class/order/family proximity, locality-count imbalance, and exact source/target identity.

The legacy v0.3 state remains `NOT_EVALUABLE_B`. TTF-Q may characterize any response-blind B relation geometry that can be reconstructed, but such characterization is method development and must never be reported as the missing v0.3 confirmatory test.

### C — historical climate displacement similarity

The focal relation is `z_R_hist`, with `z_R_current` added as the first control before the same geography/lineage/locality/source/target structure.

The legacy v0.3 state remains `NOT_EVALUABLE`. The bounded occurrence run left unresolved transport errors even though hundreds of species had usable occurrence geometry. TTF-Q may use response-blind external-data geometry in a new method-development analysis, but it cannot rewrite the terminal v0.3 receipt.

C is especially informative because its non-redundancy axis has direct ecological meaning: present-day niche similarity can arise from similar or different late-Quaternary climatic trajectories.

## Interpretation boundary

TTF-Q v0.1 is intentionally not another predictor search.

It does not:

- open nucleotide identity, pairwise genetic distances, `T_st`, `beta_R`, or `beta_hist`;
- replace B or C in relational v0.3;
- recycle the v0.3 alpha allocation;
- call transport failure a biological null;
- call a synthetic power result evidence that an ecological effect exists;
- convert its continuous diagnostics back into a post-hoc binary success gate.

The closed v0.3 audit remains immutable. TTF-Q is a fresh methods program built from the response-blind part of the machinery.

## Implementation

Machine-readable contract:

`docs/supporting/ttf_q_v0.1.json`

Reusable implementation:

`src/ttf/relational_qualification.py`

Case-study runner:

`scripts/run_ttf_q_characterization.py`

Example:

```bash
python scripts/run_ttf_q_characterization.py \
  --design path/to/relational-historical-opportunity.npz \
  --kind historical \
  --panel development \
  --output benchmarks/development/ttf_q_historical_characterization_v0.1.json
```

If response-blind bootstrap relation replicates are available, pass a NumPy matrix with replicates in rows and fixed dyads in columns using `--relation-replicates`.

## Next empirical step

The priority case study is C. Reconstruct the historical/current relation geometry from the already response-blind external assets under a new TTF-Q provenance chain, then quantify:

1. `R_hist` repeatability under occurrence resampling;
2. the fraction of `R_hist` surviving `R_current` and all frozen controls;
3. whether unique historical information is broadly distributed or concentrated in a few endpoints;
4. the effect-size × heterogeneity detectability surface.

If the historical relation contains substantial, broadly distributed information after conditioning on current climate, that is a positive ecological result about historical climatic information even before any genetic-response analysis is considered.
