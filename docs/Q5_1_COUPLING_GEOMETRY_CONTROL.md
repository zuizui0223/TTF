# Q5.1 — training-only geometry control for TTF-C

Status: **prospective successor frozen before Q5.1 qualification outcomes**.

## Why Q5.1 exists

The frozen Q5 result remains authoritative and failed jointly:

- TTF-M passed Q5;
- TTF-C failed Q5;
- joint Q5 failed.

The subsequent response-blind factor-isolation diagnostic did not select a new method. It showed a narrower mechanism pattern: TTF-C false-positive inflation was most strongly associated with clustered sampling geometry, especially clustering shift, while shared-mismatch power loss was largest under the combined heterogeneous geometry. Sample-size heterogeneity alone and response-blind missingness alone produced little change from balanced geometry.

Those observations are diagnostic only. They do not alter Q5 and do not authorize empirical use.

## Why this candidate is not invented from scratch

Core TTF had already encountered the same broad nuisance pathway during v0.2-v0.3 development: irregular sampling produces spatially structured kNN edge lengths, and edge span is associated with turnover even when biological transition locations are species-private.

The frozen v0.3 ablation compared three response-blind corrections. The least invasive candidate that met the development thresholds was **training-only edge-length orthogonalization**. It reduced private false transfer while retaining more shared-signal power than evaluation-side or two-sided conditioning.

Q5.1 therefore imports that already-defined correction rather than searching a new correction family after seeing Q5.

## Estimator change

TTF-M is unchanged.

For TTF-C only, let `u^C_se` be the ordinary within-system coupling-turnover rank and `l_se` the corresponding kNN edge-length rank. In each training system, Q5.1 replaces `u^C` by the centered rank residual after linear projection on `l`, then rescales the residual to the population SD of an equally spaced rank01 vector of the same edge count.

The transferable coupling field is fit to these training residuals with prior mean zero.

Nothing changes on the held-out side:

- same coordinates;
- same density-scaled kNN graph;
- same raw coupling-turnover rank;
- same edge-integrated field exposure;
- same raw within-system Spearman score;
- same equal-system macro-average;
- same centered held-out-system bootstrap.

Thus Q5.1 removes only a training-side geometry path. It does not partial out held-out outcomes or select a correction using held-out performance.

## Frozen formal gate

The exact rule is `docs/supporting/ttfm_q5_1_coupling_geometry_rule.json`.

The qualification repeats the seven Q5 roles on fresh worlds under master seed `20260917`, 500 worlds per cell and 1,999 bootstrap resamples per world. Numerical gates are unchanged:

- every declared type-I Wilson 95% upper bound must be <= 0.10;
- every declared power Wilson 95% lower bound must be >= 0.80.

The matched and shifted heterogeneous geometry families, graph fraction 0.15, bandwidth 0.2, prior strength 0.25 and segment quadrature are unchanged from Q5.

TTF-M is recomputed only as an invariance audit and must retain its declared roles. Raw uncorrected TTF-C is retained in the batch output for audit but cannot determine the Q5.1 decision.

## Interpretation firewall

A Q5.1 PASS would license only the statement that this fixed TTF-C successor controls the declared private-transition errors and detects the declared shared transitions on the finite frozen Q5.1 geometry family.

It would not:

- rescue or reinterpret the frozen Q5 failure;
- validate flower-colour sharedness;
- establish universal geometry robustness;
- establish a causal effect of clustering;
- establish fitness loss or mechanism breakdown;
- bypass an application-specific observation-support gate.

A Q5.1 FAIL is also terminal for this candidate. No threshold, graph fraction, bandwidth, residualization side, or qualification cell may be changed after outcomes are opened. Any successor would require a newly named method, fresh rule, and fresh worlds.
