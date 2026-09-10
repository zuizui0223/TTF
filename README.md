# TTF — Transferable Turnover Fields

TTF is a standalone, trait-agnostic methods repository for testing whether **within-species trait-transition structure learned from some species predicts transition structure in entirely unseen species**.

This repository is intentionally separated from `fcp`. `fcp` remains an empirical flower-colour ecology project; TTF is method development and qualification.

## Core claim

TTF does **not** define cross-species sharedness as hotspot concentration in one map.

> **Sharedness = out-of-species transferability.**

For species \(s\), TTF builds a graph only among records of that species, measures edge-wise trait dissimilarity, and rank-standardizes it within species. A boundary field is learned from training species, then scored on disjoint evaluation species at the edge level:

\[
C_s = \operatorname{Spearman}(\hat b_{se}, u_{se}), \qquad
T = |S_{\rm eval}|^{-1}\sum_s C_s.
\]

The method has been developed under prospective synthetic, fixed-geometry, nuisance-profiled, and observation-support gates. Earlier failed versions remain part of the frozen qualification history rather than being reclassified after later improvements.

## Qualification status

### Idealized synthetic qualification

The v0.2 held-out-species inference passed its high-precision idealized synthetic gate with controlled private-world rejection and high shared-world power. That result did **not** license arbitrary empirical deployment.

### Empirical-geometry development

Subsequent real-geometry stress tests exposed two distinct problems: private spatial transitions can be confounded with shared structure by sampling geometry, and fixed-k graphs can become too local as sampling density increases. The v0.10 geometry-only audit showed that increasing records from 20 to 100 under fixed `k=3` reduced normalized edge scale to about **0.319** and shared-transition observability to about **0.335**, whereas keeping `k/n ≈ 0.15` retained both at about **0.867**.

### v0.11 fresh confirmatory method qualification — PASS

A fresh, species-disjoint plant panel was prospectively frozen at **250 species × 100 records = 25,000 records**, with a single predeclared density-scaled graph choice `k=15`, a 125/125 train/evaluation split, 500 observed worlds and 1,999 private-reference worlds.

The confirmatory result passed both frozen gates:

- private rejection at amplitudes 0.5 / 1 / 2 / 3: **0.046 / 0.070 / 0.020 / 0.050**;
- maximum Wilson 95% upper bound: **0.09580 ≤ 0.10**;
- shared amplitude-2 power: **0.998**;
- Wilson 95% lower bound: **0.98876 ≥ 0.80**.

The exact dense operator would require roughly 103.6 GiB on that geometry, so TTF uses a chunked exact implementation that was numerically checked against the dense scorer to `1e-12`; this changes execution, not the estimand.

Frozen evidence is in `results/v11_fresh_density_scaled_qualification_v0.1.json` and the associated frozen geometry ledgers.

## v0.12 actual observation-support gate — NOT EVALUABLE

Passing v0.11 did not automatically authorize empirical flower-colour inference. The same frozen 25,000 photos were therefore sent through an independently frozen, location-blind measurement pipeline before TTF was permitted to read colour vectors, pairwise colour distances, edge turnover, or the empirical transfer statistic.

The prospective requirement was that **all 250 species retain at least 40 evaluable photographs**. The terminal result was:

- species meeting the minimum: **189 / 250**;
- median evaluable photographs per species: **51**;
- minimum observed evaluable photographs for a species: **15**;
- retained actual geometry: **0 species / 0 records**;
- `colour_vector_values_read_by_ttf = false`;
- `pairwise_colour_distances_computed = false`;
- `ttf_empirical_statistic_computed = false`.

Therefore v0.12 is a **measurement-support / observational-admissibility failure**, not an ecological null and not evidence that flower-colour sharedness is absent. The frozen rule requires stopping without species replacement, extra photographs, target relaxation, or outcome-driven tuning.

The terminal ledger is `benchmarks/frozen/v12_actual_geometry_source.json`, frozen by commit `3e4b221f2951aaf940d4695288254e27404fe2fb`.

## Paired-state extension: TTF-M and TTF-C

The separate branch `method/ttfm-mismatch-decoupling` extends the same held-out-system transfer idea to paired ecological states without altering the frozen v0.1-v0.12 history.

- **TTF-M:** transfer turnover of a predeclared mismatch magnitude `M_s(x)=m(A_s(x),B_s(x))`.
- **TTF-C:** transfer turnover of a predeclared relation state `R_s(x)=r(A_s(x),B_s(x))`.

The direction-free idealized qualification passed for both estimands. Under a shared relation rotation with exactly constant mismatch magnitude, TTF-M rejected `0/500` while TTF-C rejected `500/500`, showing that mismatch-magnitude and relation-transition fronts are distinct implemented targets.

### Q4 directionality

A first bounded-center Q4 v0.1 protocol is preserved as an **immutable formal FAIL**. The exact private front locations were independent but restricted to a common central zone, which created transferable turnover-intensity structure. That failure motivated a new protocol rather than threshold retuning.

Q4.1 uses system-private phase anchors uniformly over a periodic domain while retaining a separate outcome-independent local orientation coordinate. Q4.1 **passed** its frozen gate: shared breakdown and recoupling were labelled in `500/500` worlds, neutral rotation and cross-system sign conflict produced no directional label, and the strict-private cell had relation/breakdown rejection `24/500 = 0.048` with Wilson upper `0.07043`.

`breakdown` remains a bounded statistical label: transferable increase in the predeclared mismatch coordinate across a predeclared signed orientation. It is not by itself causality, fitness loss, or mechanism failure.

### Q5 heterogeneous geometry stress — ACTIVE FORMAL GATE

Q5 is now authorized and running on a finite predeclared family of response-blind geometries. It varies nominal records per system (`30/45/60/90/120`), clustering, domain truncation, `0/10/25%` missingness, and matched versus shifted train/evaluation geometry. The graph degree is selected separately for each realized system by the frozen rule `k_s=max(2, round(0.15*n_s))` using post-missingness effective `n_s`.

Seven mandatory cells test component-only nulls, private mismatch fronts, private relation rotations, shared mismatch/relation power, and lower-density/more-clustered evaluation shifts. Formal outcomes remain unopened until all batches aggregate. Q5 PASS would qualify only this finite geometry family and would still require an independent application-specific observation-support gate before any empirical use.

## Why the observation-support distinction matters

TTF separates three claims that must not be collapsed:

```text
method validity on declared worlds
        !=
deployment-geometry detectability
        !=
actual observational admissibility
```

A method may be statistically valid and powerful on a prospectively qualified geometry while the realized measurement process fails to preserve enough support to instantiate that geometry empirically. In that case the correct endpoint is **not evaluable / abstain**, not a biological negative.

## Install and test

```bash
python -m pip install -e ".[test]"
pytest -q
```

## Claim ceiling

TTF currently supports bounded method statements only. The original flower-colour route remains non-evaluable before any empirical TTF statistic was opened. The TTF-M/TTF-C extension is qualified on its declared balanced idealized family and Q4.1 directionality family; Q5 heterogeneous geometry qualification is still active. No empirical flower-colour, pollinator, island, fitness, causal mismatch, or universal robustness conclusion is licensed.
