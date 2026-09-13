# Rank-symmetry obstruction for residual TTF-C conditional-null bias

## Scope

This note is a deductive consequence of the merged phase-propensity diagnosis and the exact rank-factorization result. It does **not** define a new estimator, repair candidate, mechanism-development run, or qualification route.

The frozen phase-propensity result remains immutable: subtracting the exact geometry-only phase-expected unscaled training response removed most of the positive heterogeneous conditional-null mean, but a smaller positive residual remained. Same-world crossing-propensity, angular-distance, residual-offset, normalization, leverage, and support-concentration tuning remain closed.

## 1. The centered arm has zero raw predicted mean

In the frozen phase-propensity centered arm, each training edge response has its exact conditional phase mean subtracted. Therefore, conditional on fixed geometry `G`,

`E[T | G] = 0`

edge by edge for the centered training response `T`.

The frozen transfer operator is geometry-only and linear in `T`. In `score_chunked_batch`, with `prior_mean = 0`, the held-out predicted edge-value vector is

`P = A_G T`,

where `A_G` is fixed once graph geometry, bandwidth, prior strength, and quadrature are fixed. Hence

`E[P | G] = A_G E[T | G] = 0`.

So the residual positive conditional-null Spearman expectation does not require any remaining raw prediction mean.

## 2. Rank direction is odd

Let `u(v)` be the standardized centered average-rank direction used in the exact rank-factorization note, with constant vectors mapped to zero.

For every finite vector `v`, including vectors with ties,

`u(-v) = -u(v)`.

Average ranks reverse under global sign reversal, centered ranks change sign, and their Euclidean norm is unchanged. Constant vectors remain the zero vector.

## 3. Central symmetry forces zero expected rank direction

If the conditional prediction law is centrally symmetric,

`P | G  =d=  -P | G`,

then oddness gives immediately

`E[u(P) | G] = 0`.

This conclusion is unaffected by unequal marginal variances, covariance structure, or heteroskedasticity across held-out edges. In particular, any zero-centered multivariate Gaussian or other centrally symmetric elliptical prediction law has zero expected rank direction even when its covariance is highly nonuniform.

Combined with the merged exact factorization

`E[rho_S(P,Y) | G] = E[u(P) | G]^T E[u(Y) | G]`,

central symmetry of `P | G` would force

`E[rho_S(P,Y) | G] = 0`.

## 4. Consequence of the frozen positive residual

The phase-propensity centered arm nevertheless retained a prospectively frozen positive conditional-null mean on both heterogeneous cells:

- matched geometry: `0.0045178`, 95% CI `[0.0033935, 0.0055877]`;
- shifted geometry: `0.0065181`, 95% CI `[0.0055625, 0.0073991]`.

Therefore the residual cannot be explained by a **purely variance/covariance-based, centrally symmetric zero-mean prediction law**. The conditional prediction distribution must retain rank-relevant distributional asymmetry: equivalently, its realizations cannot cancel under global sign reversal in the way a centrally symmetric law would.

This is a stronger localization than “Spearman is nonlinear.” The residual is constrained to a higher-order ordering asymmetry of the zero-mean predicted field, aligned with the held-out expected rank field.

## 5. Zero raw mean is not enough

The converse fails: `E[P | G] = 0` does not imply `E[u(P) | G] = 0`. A finite asymmetric support can have exactly zero coordinate-wise raw mean while its expected rank direction is nonzero. The accompanying tests include such a construction.

Thus the already-frozen analytic mean subtraction can remove the first moment yet leave a nonzero rank expectation through asymmetric higher-order structure.

## What this does not establish

This note does not identify which geometric feature creates the remaining asymmetry, and it does not authorize skewness correction, pairwise-win correction, rank centering, variance normalization, or any other repair on the phase-propensity worlds.

Any empirical or simulation-based attempt to explain the residual ordering asymmetry is a **new mechanism family** and requires:

- a new name;
- a prospectively frozen prediction;
- fresh geometry worlds and phase seeds;
- no tuning or candidate selection on any already-opened diagnostic worlds.

Heterogeneous-geometry TTF-C remains unqualified. `candidate_selected = null` and `formal_pass_fail = null` remain unchanged.
