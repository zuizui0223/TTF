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
- generic held-out predictor-space competition with paired incremental skill;
- **Gate-I fixed-geometry semi-synthetic calibration**, which retains empirical species identities, coordinates, record counts and optional dependence blocks while replacing all trait values with synthetic transitions;
- SHA-256 fingerprinting of the exact empirical sampling frame used for Gate-I qualification.

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
  --replicates 100 \
  --resamples 1999 \
  --inference heldout_species_bootstrap \
  --output calibration.json
```

## Gate I — qualify the intended empirical geometry

Gate I keeps the empirical sampling frame fixed and simulates only the trait. The CSV may contain empirical trait columns, but `run_geometry_calibration.py` never reads them. Coordinates should already be in the coordinate system intended for TTF graph construction and bandwidth selection.

```bash
python scripts/run_geometry_calibration.py \
  --input benchmark_geometry.csv \
  --species-column species \
  --coordinate-columns x,y \
  --bandwidth 25 \
  --shared-fractions 0,1 \
  --amplitudes 0.5,1,2,3 \
  --replicates 500 \
  --resamples 1999 \
  --output results/gate_i_benchmark.json
```

The output freezes the geometry fingerprint, per-species record counts, the prospective train/evaluation split, all simulation/inference settings, point estimates, and Wilson-bound qualification. An independent non-flower empirical sampling geometry should be preferred for the first claim-bearing Gate-I run.

**Current status:** the Gate-I engine is implemented, but no independent empirical geometry has yet been frozen and passed at high precision. The method therefore remains qualified only for its idealized synthetic scope until that run is completed.

Empirical deployment should not precede qualification on both synthetic worlds and the intended sampling geometry.

## Claim ceiling

Passing TTF can support a statement of the form:

> A transition structure learned from one set of species predicts where strong within-species trait transitions occur in unseen species, within the measurable sampling frame.

It does not by itself prove a universal causal barrier, eliminate observation bias, or identify a mechanism. Predictor-space attribution, IBD residualization, and observation-bias controls are separate qualification layers.
