# Independent historical-signal contrast: oracle design benchmark

Date: 2026-11-14. This is a separate, deterministic methods benchmark, not v2 approval of the closed Haldane test. No empirical rates, no new simulated worlds, no revised previous PASS/FAIL.

## Question and assumptions

Distinguish prediction gain of a fixed regularized learner from the coefficient of one history term conditional on the stated nuisance span. Consider y=Z beta + gamma h + epsilon, epsilon ~ N(0, Sigma), where the FULL covariance Sigma is specified correctly in advance. This assumption is an oracle benchmark, NOT available for actual PROCEED trait rates. Neither normality nor real variance components have been verified.

Use the original source-only arrays only if their original byte hashes match. Keep the 29/27 split. The hypothetical test uses only the 74 heldout source-pair rows; training rows fix the scale of h (equal total weight per training species). It is a heldout conditional-association design, NOT an out-of-species prediction test. A common gamma is assumed, not established.

## Options and chosen scope

1. Refit/tune ridge or enlarge old Monte Carlo calibration: rejected; the old family is closed.
2. Study-species component resampling: potentially useful future work, but requires separate calibration and assumptions; not implemented here.
3. Analytic GLS nuisance-orthogonal contrast with known covariance: chosen as an auditable ideal benchmark. It exposes which independence/calibration assumptions must be supplied before any empirical test, without changing the old family.

## Fixed math and reporting

Cholesky Sigma=L L^T. Whiten Z and h. Project whitened h off the column space of whitened Z using scaled-column SVD, so duplicated nuisance columns do not change the contrast. Let r be this residual and I=r^T r. If I is numerically zero, report nonidentifiable. Otherwise u=L^(-T)r/sqrt(I), u^T Z=0, u^T Sigma u=1, and u^T h=sqrt(I). Under the stated oracle null u^T y ~ N(0,1), irrespective of nuisance beta.

No y is read or generated. Integrate normal tails directly for the fixed alpha=1/60 and gamma grid {0,.25,.5,1}; gamma is in hypothetical response units per original train SD of h. Report the gamma required for 80% analytic power as a DESIGN BENCHMARK, not an estimated effect. The primary known covariance retains the original diagnostic .25*(species-match + study-match + identity). A prespecified diagonal-only assumed covariance is compared against that same true covariance to illustrate sensitivity, without choosing between models based on results. Neither is a biological estimate.

Hypothetical association calibration is model-specific; mean-model errors, estimated covariance, trait-rate skew, phylogeny, ancestral occupancy, and heterogeneous slopes remain unresolved. The known-covariance benchmark cannot be used to authorise empirical opening or claim adaptation/plasticity/speciation.

## Implementation plan

- [ ] Write deterministic tests for nuisance annihilation, unit variance, duplicate history rejection, duplicate-control invariance, study covariance, bad inputs, normal-tail size and power, zero random draws.
- [ ] Implement pure source-design contrast/power functions without a response-input or RNG interface.
- [ ] Run local tests including invariance and one-ULP strict source guards. Independent review of algebra and input handling.
- [ ] Run the source-only benchmark on GitHub Actions after original code, native CSV, coordinates, climate, X and h hashes pass. Original helper files remain unchanged. Read source metadata only; no real means/SD/rates parsed.
- [ ] Export ONLY aggregate information, analytic size/power and source hashes. No contrasts, original arrays, pointwise coordinates or labels are uploaded.
- [ ] Verify the original terminal record and empirical STOP test unchanged; save result with clear oracle conditions.

References: statsmodels GLS official documentation (Sigma and whitening); SciPy norm official documentation (normal cdf/sf/isf). Original TTF diagnostic commit: 53f391a5873a31db0e83cb952bd801a10211b932; original source model commit: a4664cc1fccf7a16fed1b85f239dc0f123389a87.
