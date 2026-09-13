# Exact rank-factorization of residual TTF-C conditional-null bias

## Scope

This note does **not** define a new estimator, mechanism-development family, or qualification route. It is an algebraic localization result for the residual conditional-null bias that remains after the frozen edge-phase propensity diagnosis.

The frozen phase-propensity result is immutable. Its analytic training-response subtraction removed about 78% of the matched and 71% of the shifted positive conditional-null mean, while a smaller positive residual remained. No further crossing-propensity, angular-distance, residual-offset, normalization, leverage, or support-concentration tuning is licensed on those worlds.

## Spearman as a rank-direction inner product

For a nonconstant vector `v`, let

`u(v) = (R(v) - mean(R(v))) / ||R(v) - mean(R(v))||`,

where `R(v)` is the average-rank vector with deterministic tie handling. For a constant vector define `u(v)=0`, matching the repository's `spearman_rho` contract.

Then exactly

`rho_S(x,y) = u(x)^T u(y)`.

No approximation is involved.

## Conditional independence factorization

Fix all graph and sampling geometry `G` for one held-out system. After the phase-propensity training-response mean has been subtracted, let

- `P` be the predicted edge-value vector generated from the independently drawn training-system private phases;
- `Y` be the held-out coupling-turnover vector generated from that held-out system's private phase.

Under the frozen private-phase design, the training phases and held-out phase are independent conditional on `G`. Therefore `u(P)` and `u(Y)` are conditionally independent. Hence

`E[rho_S(P,Y) | G]`

`= E[u(P)^T u(Y) | G]`

`= E[u(P) | G]^T E[u(Y) | G]`.

This identity remains valid with ties and with constant realizations because the zero-vector convention is part of `u`.

## Consequence for the frozen residual

Once the geometry-only mean training response has been removed, a positive conditional-null Spearman expectation cannot require a remaining nonzero raw training-response mean. It must live in the **rank layer**:

1. geometry makes some evaluation edges systematically occupy different expected normalized rank positions under random private phase;
2. geometry also makes the rank ordering of the zero-mean predicted field non-exchangeable across evaluation edges;
3. the conditional-null expectation is exactly the alignment of those two expected normalized rank fields.

Thus the remaining phase-propensity residual is mathematically localized to **geometry-conditioned expected-rank alignment**. This is stronger than saying merely that Spearman is nonlinear: it identifies the exact object whose inner product equals the residual expectation.

## What this does not establish

The factorization is an identity, not a prospectively tested biological or statistical mechanism. It does not establish why the expected predicted-rank field is nonzero, nor does it authorize a rank-centering correction.

Any attempt to explain or remove that expected-rank field is a **new mechanism family** and, under the repository governance, requires:

- a new name;
- a new prospectively frozen prediction;
- fresh geometry worlds and phase seeds;
- no reuse of the phase-propensity worlds for tuning or candidate selection.

Heterogeneous-geometry TTF-C therefore remains unqualified, and `candidate_selected = null` / `formal_pass_fail = null` remain the governing endpoint.
