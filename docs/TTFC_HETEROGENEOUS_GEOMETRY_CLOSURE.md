# TTF-C heterogeneous-geometry closure ledger

## Status

**Terminal methodological closure for the current TTF-C heterogeneous-geometry estimator family.**

This document consolidates the frozen Q5/Q5.1/Q5.2 and subsequent diagnostic chain. It does not reopen any failed gate, select a repair, or authorize reuse of diagnostic worlds.

Current endpoint:

- heterogeneous-geometry TTF-M: qualified under Q5;
- heterogeneous-geometry TTF-C: **not qualified**;
- `candidate_selected = null`;
- `formal_pass_fail = null` for diagnostic families;
- no empirical flower-colour conclusion is licensed by these synthetic diagnostics.

## 1. What failed

Q5 established that heterogeneous observation geometry breaks the joint TTF-M/TTF-C qualification even though TTF-M itself remains robust. Q5.1 showed that training-side edge-length residualization is only a partial control for TTF-C:

- private false transfer remains too high;
- shifted shared-mismatch power remains inadequate;
- pure shared relation remains easy;
- therefore the failure is not a generic inability to detect relational turnover.

The scientific problem is narrower: under heterogeneous geometry, the frozen TTF-C estimator acquires a positive conditional-null expectation for private relation phases.

## 2. Mechanisms prospectively tested and closed

### Support overlap — not supported

Simple kernel opportunity / support overlap did not explain the residual heterogeneous TTF-C failure. No support-threshold or opportunity tuning is licensed on those worlds.

Frozen result: `results/ttfc_support_overlap_diagnostic_v0.1.json`.

### Training-system leverage — not supported as primary

Deletion sensitivity is real and can be stronger for TTF-C, but leverage/concentration did not explain private-C false transfer. One estimand-specific shared-mismatch effect was retained as a bounded secondary observation only.

Frozen result: `results/ttfc_training_leverage_diagnostic_v0.1.json`.

Leverage and support-concentration tuning are closed.

### Conditional-null centering — supported

With geometry fixed and private relation phases redrawn prospectively, heterogeneous Q5 geometry produced a reproducible positive phase-averaged TTF-C statistic while the homogeneous control was near zero.

This established that the residual false transfer is not an accidental alignment of one realized private phase. It is a geometry-conditioned null expectation of the frozen estimator.

Frozen result: `results/ttfc_conditional_null_centering_v0.1.json`.

### Final residual-SD normalization — contributor, not primary

Removing the final residual-SD rescaling reduced the positive null bias but did not recenter it.

Matched heterogeneous geometry:

- normalized mean: `0.0231777`;
- unscaled mean: `0.0211155`;
- normalization contribution: `0.0020622`.

Shifted heterogeneous geometry:

- normalized mean: `0.0238779`;
- unscaled mean: `0.0223741`;
- normalization contribution: `0.0015037`.

Thus normalization amplifies the problem but is not necessary for it. Normalization floors, clipping, shrinkage, and alternative scaling are closed on those worlds.

Frozen result: `results/ttfc_normalization_bias_diagnostic_v0.1.json`.

### Edge-phase crossing propensity — major contributor, not complete mechanism

For a Q5 private relation boundary with uniform phase, an edge with shortest endpoint angular separation `delta_e` crosses with exact probability

`p_e = delta_e / pi`.

The exact geometry-only phase expectation of the unscaled Q5.1 training response was derived analytically and subtracted without fitted parameters.

Matched heterogeneous geometry:

- baseline unscaled conditional mean: `0.0209832`;
- baseline minus centered: `0.0164654`;
- centered residual: `0.0045178`, 95% CI `[0.0033935, 0.0055877]`.

The analytic mean accounts for about **78.5%** of the baseline mean, but the residual remains positive.

Shifted heterogeneous geometry:

- baseline unscaled conditional mean: `0.0223623`;
- baseline minus centered: `0.0158442`;
- centered residual: `0.0065181`, 95% CI `[0.0055625, 0.0073991]`.

The analytic mean accounts for about **70.9%** of the baseline mean, but the residual remains positive.

Therefore the prospectively frozen necessary-mechanism conjunction failed: edge-phase propensity is a major upstream contributor, not the complete mechanism. Crossing-propensity centering, angular-distance transforms, residual offsets, and related threshold tuning are closed on those worlds.

Frozen result: `results/ttfc_phase_propensity_diagnostic_v0.1.json`.

## 3. Exact mathematical localization of the residual

The repository's Spearman score can be written exactly as an inner product of standardized centered rank directions:

`rho_S(P,Y) = u(P)^T u(Y)`.

Because training private phases and a held-out private phase are independent conditional on fixed geometry,

`E[rho_S(P,Y) | G] = E[u(P) | G]^T E[u(Y) | G]`.

Therefore the residual after analytic mean subtraction is not an unidentified scalar nuisance. It is exactly an alignment between two geometry-conditioned expected rank fields.

See `docs/TTFC_RANK_FACTORIZATION.md`.

## 4. Symmetry obstruction: what the residual cannot be

In the phase-propensity centered arm, the training response has exact conditional phase mean zero edge by edge. With `prior_mean = 0`, the frozen transfer estimator is a geometry-fixed linear map, so the held-out predicted edge-value vector also has

`E[P | G] = 0`.

The standardized rank direction is odd, including under ties:

`u(-v) = -u(v)`.

Hence any centrally symmetric conditional prediction law

`P | G =d= -P | G`

must satisfy

`E[u(P) | G] = 0`,

which would force zero expected Spearman under the exact factorization above.

The frozen positive residual therefore rules out a **purely variance/covariance-based centrally symmetric zero-mean account**. Unequal variance or covariance alone, within a centrally symmetric family such as a zero-centered multivariate Gaussian or other elliptical law, cannot generate the observed residual expectation.

What remains is rank-relevant **distributional ordering asymmetry** of the zero-mean predicted field. Zero first moment is not sufficient: an asymmetric finite support can have exact coordinate-wise mean zero while retaining a nonzero expected rank direction.

See `docs/TTFC_RANK_SYMMETRY_OBSTRUCTION.md`.

## 5. Scientific interpretation

The heterogeneous-geometry result is now asymmetric between the two TTF estimands:

- **TTF-M** can transfer mismatch magnitude under the frozen heterogeneous-geometry stress test.
- **TTF-C** is more fragile because relational turnover interacts with observation geometry to create estimator-induced null structure even when private relation phases are independent across systems.

This is itself a useful methodological conclusion. It means magnitude sharedness and relational/coupling sharedness should not be treated as equally identifiable under heterogeneous sampling geometry.

The appropriate paper-level interpretation is therefore not “TTF-C almost works and needs one more correction.” It is:

> Under heterogeneous geometry, the current TTF-C estimator has a demonstrable geometry-conditioned identifiability failure. Most of that bias is analytically attributable to edge-phase crossing propensity, normalization adds a smaller amplification, and the remaining component is constrained to non-centrally-symmetric rank-ordering structure. No prospectively qualified correction has been selected.

## 6. Stop rule

For the current manuscript and estimator family, **stop TTF-C mechanism hunting here** unless a new theoretical result creates a qualitatively different, prospectively falsifiable mechanism with clear information value.

Specifically, do not open new same-purpose families merely to reduce the remaining `~0.0045–0.0065` residual by trying:

- rank centering;
- skewness correction;
- pairwise-win correction;
- alternate angular transforms;
- variance stabilization;
- residual offsets;
- support/leverage reweighting;
- normalization variants.

Any future simulation-based successor must have all of the following before outcomes are generated:

1. a new mechanism name;
2. a mechanistic prediction not equivalent to an already-closed family;
3. a prospectively frozen pass/fail rule;
4. fresh geometry worlds and phase seeds;
5. no reuse of opened worlds for tuning or candidate selection.

Absent that, the correct endpoint is the present bounded failure, not another rescue attempt.

## 7. Manuscript-facing claim boundary

Allowed:

- report TTF-M heterogeneous-geometry robustness under its frozen gate;
- report that TTF-C fails heterogeneous-geometry qualification;
- report the supported geometry-conditioned null bias;
- report normalization as a small amplifier;
- report edge-phase propensity as a major but incomplete contributor;
- use the exact rank factorization and symmetry obstruction to localize the unresolved remainder;
- present this as an identifiability/estimator limitation.

Not allowed:

- claim a qualified heterogeneous-geometry TTF-C estimator;
- promote the propensity-centered ablation to a selected estimator;
- claim the remaining bias is biologically causal;
- infer empirical flower-colour sharedness or no-sharedness from these synthetic diagnostics;
- tune a new correction on any diagnostic worlds already opened.
