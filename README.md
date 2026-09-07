# TTF — Transferable Turnover Fields

TTF is a standalone methods repository for testing whether **within-species trait-transition structure learned from some species predicts transition structure in entirely unseen species**.

This repository is intentionally separated from `fcp`. `fcp` remains an empirical flower-colour ecology project; TTF is trait-agnostic method development.

## Core claim

TTF does **not** define cross-species sharedness as hotspot concentration in one map.

> **Sharedness = out-of-species transferability.**

For species \(s\), TTF builds a graph only among records of that species, measures edge-wise trait dissimilarity, and rank-standardizes it within species. A boundary field is learned from training species, then scored on held-out species at the **edge level**:

\[
C_s = \operatorname{Spearman}(\hat b_{se}, u_{se}), \qquad
T = |S_{\rm eval}|^{-1}\sum_s C_s.
\]

Here \(\hat b_{se}\) is boundary exposure learned without species \(s\), integrated along edge \(e\), and \(u_{se}\) is the within-species turnover rank.

## Why this repository exists

The RGFCA experience exposed two quantities that must not be conflated:

1. **within-species spatial amplitude** — a species can have a very strong geographic transition;
2. **cross-species sharedness** — independent species put transitions in the same places.

A method can falsely call (1) evidence for (2). TTF therefore requires a calibration arm where every species has strong spatial structure but **shared fraction = 0**.

## Implemented in the first methods milestone

- within-species-only kNN graphs;
- generic scalar/vector trait dissimilarity;
- within-species rank standardization;
- species-equal opportunity-corrected kernel boundary fields with shrinkage;
- balanced repeated inclusion schedules;
- descriptive recurrence probability (kept separate from inference);
- species-disjoint edge-level transfer statistic;
- full within-species trait-location permutation null that refits the training field every replicate;
- optional block-restricted permutation for exchangeability;
- geometry-only kernel precomputation that is numerically equivalent to direct scoring;
- synthetic shared-fraction × amplitude calibration worlds;
- generic held-out predictor-space competition with paired incremental skill.

See [`docs/METHOD_SPEC.md`](docs/METHOD_SPEC.md) and [`docs/QUALIFICATION_PROTOCOL.md`](docs/QUALIFICATION_PROTOCOL.md).

## Install and test

```bash
python -m pip install -e ".[test]"
pytest -q
```

## Run a calibration grid

```bash
python scripts/run_calibration.py \
  --shared-fractions 0,0.25,0.5,1 \
  --amplitudes 0,1,2,3 \
  --replicates 25 \
  --permutations 199 \
  --output calibration.json
```

The calibration command is deliberately explicit. Empirical deployment should not precede qualification on both idealized synthetic worlds and the actual sampling geometry.

## Claim ceiling

Passing TTF can support a statement of the form:

> A transition structure learned from one set of species predicts where strong within-species trait transitions occur in disjoint species, within the measurable sampling frame.

It does not by itself prove a universal causal barrier, nor does it eliminate observation bias. Predictor-space attribution and negative-control sampling analyses are separate layers.
