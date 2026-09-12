# TTF v0.3 edge-geometry development record

## Status

v0.3 is a **successor-development line**, not a reinterpretation of v0.2.
The v0.2 Gate-I failures remain frozen and retain their original decisions.

The selected v0.3 candidate is:

> **train_orthogonalized_only** — within each training species, project the
> within-species turnover rank off the within-species kNN edge-length rank before
> fitting the boundary field. Held-out traits and the held-out raw Spearman score
> are unchanged.

This candidate is not qualified until it passes the prospectively frozen fresh
external-geometry qualification described below.

## Why v0.2 Gate I failed

### Gate I-A: independent animal geometry

The frozen 16-species external animal panel used 8 training and 8 held-out
species. At 500 worlds per cell and 1,999 held-out-species bootstrap resamples:

- shared=0, amplitude 0.5: rejection 0.058;
- shared=0, amplitude 1: rejection 0.046;
- shared=0, amplitude 2: rejection 0.072;
- shared=0, amplitude 3: rejection **0.110**, Wilson 95% upper **0.1405**;
- shared=1, amplitude 2: power **0.392**, Wilson 95% lower **0.3502**.

Thus I-A failed both the type-I and power gates.

### Gate I-B: intended RGFCA geometry

The intended geometry was reconstructed without empirical colour from the exact
pre-G1/G3 FCP state `e6ad7aa76261906b6f757ffc324ea6cd062a0917`.
The frozen outer-0 frame contains 250 species x 20 photographs = 5,000 records.
The high-precision qualification found strong false transfer under private
transitions, including rejection 0.134 at shared=0/amplitude=2 and 0.210 at
shared=0/amplitude=3, while the fully shared amplitude-2 control had power 0.858.

Therefore the principal I-B problem is not lack of power. Real sampling geometry
can make strong private transition systems appear transferable.

## Edge-length nuisance mechanism

For a spatially smooth trait with a transition, the expected trait difference
between two observations increases with the spatial span of the edge even when
the transition orientation and location are private to a species. In irregular
sampling frames, kNN edge length itself is spatially structured. Similar
sampling-density geometry across species can therefore create a cross-species
predictive path:

```text
sampling geometry
      |
      +--> spatially structured kNN edge length
      |                  |
      |                  +--> larger expected turnover on long edges
      |
      +--> training boundary field
                         |
                         +--> held-out turnover association
```

This path does not require a shared biological transition mechanism.
Opportunity correction at observation locations does not by itself remove this
edge-level nuisance.

## Paired nuisance diagnostic

The failed I-B geometry was reused **only for development diagnosis**. It cannot
qualify a successor.

A candidate that controlled edge length on both training and evaluation sides
reduced private-world rejection substantially:

- shared=0, amplitude 2: approximately 0.128 -> 0.048;
- shared=0, amplitude 3: approximately 0.206 -> 0.074.

However fully shared amplitude-2 power dropped from approximately 0.862 to 0.758.
This supported the nuisance mechanism but showed that controlling both sides was
too destructive for the prospective power target.

## Frozen three-way ablation

Before the ablation outcome was inspected, three candidate corrections and a
preference order were frozen:

1. `eval_partial_only`;
2. `train_orthogonalized_only`;
3. `train_and_eval_conditioned`.

The development selection rule was the least invasive listed candidate with:

- maximum rejection <= 0.10 over shared=0, amplitude 2 and 3;
- power >= 0.80 for shared=1, amplitude 2.

Five hundred paired worlds per cell gave:

| Variant | max private rejection | shared A=2 power | development decision |
|---|---:|---:|---|
| original | 0.200 | 0.850 | fail |
| eval_partial_only | 0.090 | 0.772 | fail |
| **train_orthogonalized_only** | **0.092** | **0.824** | **select** |
| train_and_eval_conditioned | 0.080 | 0.770 | fail |

The selected candidate is therefore `train_orthogonalized_only`. No coefficient,
k, bandwidth, prior strength, or evaluation statistic was tuned to improve this
result.

## Exact v0.3 estimand change

For one training species and edge `e`, let `u_e` be the original within-species
turnover rank and `l_e` the within-species edge-length rank. v0.3 replaces the
training response by the centered rank residual

\[
  u^{\perp}_e = (u_e-\bar u) -
  \frac{\operatorname{Cov}(u,l)}{\operatorname{Var}(l)}(l_e-\bar l),
\]

then rescales the residual to the population standard deviation of an equally
spaced rank vector of the same edge count. The rescaling keeps species-level
field contribution scales comparable. Because the residual has mean zero, the
boundary-field prior mean is fixed to zero.

Everything on the evaluation side remains v0.2:

- same held-out coordinates and graph;
- same raw turnover rank target;
- same edge-integrated field exposure;
- same raw within-species Spearman score;
- same equal-species macro-average;
- same centered, studentized held-out-species bootstrap.

## Fresh external validation firewall

A 40-species North American animal panel was declared before candidate
performance was inspected. It is taxonomically disjoint from the failed I-A
panel. Acquisition freezes only occurrence geometry; no synthetic trait world is
run during geometry freeze.

The fresh geometry ledger must record:

- `confirmatory_holdout=true`;
- `candidate_performance_evaluated_at_freeze=false`;
- `synthetic_worlds_run_at_freeze=0`;
- zero taxon overlap with failed I-A;
- k=4 and a geometry-only bandwidth fixed as the median of species median kNN
  edge lengths;
- a fixed 20/20 train/evaluation split under seed 20260908.

Only after the geometry is frozen is a separate authorization file allowed to
open synthetic outcomes.

## Fresh high-precision qualification gate

The authorization fixes exactly five cells:

- shared=0, amplitude 0.5;
- shared=0, amplitude 1;
- shared=0, amplitude 2;
- shared=0, amplitude 3;
- shared=1, amplitude 2.

Each cell uses 500 worlds and 1,999 bootstrap resamples. A fresh-panel pass
requires:

- every zero-shared cell's two-sided Wilson 95% upper bound <= 0.10;
- the fully shared amplitude-2 Wilson 95% lower bound >= 0.80.

A pass qualifies only this locked v0.3 correction on the fresh external geometry.
It does not retroactively turn the v0.2 Gate-I failures into passes and does not,
by itself, qualify empirical RGFCA deployment. A subsequent intended-geometry
qualification for the locked v0.3 estimator remains a separate layer.
