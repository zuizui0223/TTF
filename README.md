# TTF — Transferable Turnover Fields

TTF is a standalone, trait-agnostic methods repository for testing whether **within-species trait-transition structure learned from some species predicts transition structure in entirely unseen species**.

This repository is intentionally separated from `fcp`. `fcp` remains an empirical flower-colour ecology project; TTF is method development.

## Core claim

TTF does **not** define cross-species sharedness as hotspot concentration in one map.

> **Sharedness = out-of-species transferability.**

For species \(s\), TTF builds a graph only among records of that species, measures edge-wise trait dissimilarity, and rank-standardizes it within species. A boundary field is learned from training species, then scored on disjoint evaluation species at the **edge level**:

\[
C_s = \operatorname{Spearman}(\hat b_{se}, u_{se}), \qquad
T = |S_{\rm eval}|^{-1}\sum_s C_s.
\]

Here \(\hat b_{se}\) is boundary exposure learned without the evaluation species, integrated along edge \(e\), and \(u_{se}\) is the within-species turnover rank.

## Primary inference in v0.2

The primary sharedness null is conditional on the frozen training field:

\[
H_0:E[C_s\mid B_{\rm train}]\le 0.
\]

TTF uses the held-out species as the resampling units in a centered, studentized species bootstrap. Coordinates, trait values, graphs, and strong private within-species spatial transitions are left untouched.

Within-species trait permutation is retained as a **diagnostic test of trait-location exchangeability**, not as the primary sharedness null.

## Why the null changed

The first prospective qualification deliberately included worlds where every species had a strong spatial transition but the transition locations were unrelated across species.

The v0.1 trait-permutation null failed that adversarial test:

- maximum zero-shared rejection = **0.14**;
- fully shared moderate-signal power = **1.00**.

The edge-level transfer statistic itself stayed centered near zero under no sharing. The failure came from a null that destroyed each species' private spatial structure.

The v0.2 held-out-species inference preserves that structure.

### High-precision idealized qualification

A frozen 500-world qualification used 1,999 held-out-species bootstrap resamples per world, 40 species and 60 records per species. For zero-shared worlds, rejection rates at amplitudes 0.5, 1, 2 and 3 were **0.040, 0.050, 0.040 and 0.054**. The worst two-sided Wilson 95% upper bound was **0.0774**, below the prospective ceiling of 0.10.

For the fully shared moderate-signal control (shared fraction 1, amplitude 2), power was **0.978** and the Wilson 95% lower bound was **0.9610**, above the prospective floor of 0.80.

Thus v0.2 **passes the high-precision type-I and power gates for its idealized synthetic scope**. This is not yet a blanket qualification for arbitrary empirical sampling geometries.

Frozen results and workflow lineage are under `results/`, including `results/qualification_heldout_species_precision_v0.2.json`.

## Implemented

- within-species-only kNN graphs;
- generic scalar/vector trait dissimilarity;
- within-species rank standardization;
- species-equal opportunity-corrected kernel boundary fields with shrinkage;
- balanced repeated inclusion schedules and descriptive recurrence;
- species-disjoint edge-level transfer statistic;
- centered, studentized held-out-species bootstrap sharedness inference;
- diagnostic within-species trait-location permutation with optional dependence blocks;
- geometry-only kernel precomputation numerically checked against direct scoring;
- synthetic shared-fraction × amplitude calibration worlds;
- confidence-bound qualification gates;
- generic held-out predictor-space competition with paired incremental skill.

See [`docs/METHOD_SPEC.md`](docs/METHOD_SPEC.md) and [`docs/QUALIFICATION_PROTOCOL.md`](docs/QUALIFICATION_PROTOCOL.md).

## PAYOFF mechanistic spatial projection

TTF can also be used as the **cross-species spatial-projection layer** for an independently specified architecture-payoff mechanism.

The current bridge assigns complementary roles:

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

A PAYOFF-derived predictor space `P` can be compared against geography/environment baselines using the existing generic predictor-space competition, provided the evaluation species' trait outcomes are not used to tune the mechanism predictor.

The Izu *C. microdonta* populations are treated as one deeply resolved anchor species, not as multiple held-out species for TTF's superpopulation claim.

See [`docs/PAYOFF_SPATIAL_PROJECTION_HANDOFF.md`](docs/PAYOFF_SPATIAL_PROJECTION_HANDOFF.md) for the feature contract, leakage rules, paired increments, directional caveat, and claim boundary.

## Run a calibration grid

```bash
python scripts/run_calibration.py \
  --shared-fractions 0,0.25,0.5,1 \
  --amplitudes 0,1,2,3 \
  --replicates 100 \
  --resamples 1999 \
  --inference heldout_species_bootstrap \
  --output calibration.json
```

Empirical deployment should not precede qualification on both synthetic worlds and the intended sampling geometry.

## Claim ceiling

Passing TTF can support a statement of the form:

> A transition structure learned from one set of species predicts where strong within-species trait transitions occur in unseen species, within the measurable sampling frame.

It does not by itself prove a universal causal barrier, eliminate observation bias, or identify a mechanism. Predictor-space attribution, IBD residualization, and observation-bias controls are separate qualification layers.
