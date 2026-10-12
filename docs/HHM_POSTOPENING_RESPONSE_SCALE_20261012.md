# HHM after outcome exposure: response magnitude versus standardized spatial pattern

Date: 2026-10-12.
Status: POST-OUTCOME ANALYTICAL NOTE, NOT A NEW CONFIRMATORY RESULT.

## What is preserved

The original one-shot 259-species / 20,362-edge standardized-response test is closed as `NOT_EVALUABLE_EMPIRICAL_RESPONSE`. Its primary statistic and Monte Carlo p-value remain null. The result file is unchanged, SHA-256 `5d994a76f7f20a778f0394f545c08a729182ce32298503c128d95c00da29f965`.

The recorded empirical run processed all 259 species and found 255 finite standardized coefficients and four undefined coefficients caused by exactly zero sampled COI variation. The four were not removed or imputed. The original outcome has already been exposed and cannot be represented as a prospective test again.

This does not imply that every possible analysis of these data is undefined. It identifies a limitation of the particular response-standardized estimand. Distinguishing estimands is essential before further work.

## Two different questions, two coefficient scales

For species s, let y_s be the endpoint-disjoint post-IBD response in unstandardized p-distance units, h_s=z(M_host), and C_s=[1,z(M_self),z(D_host_0BP)]. Define r_s=h_s-P_C h_s, where P_C projects on the controls. The already-qualified full predictor matrix has nonzero residual host information.

### A. Original standardized-pattern coefficient

For sigma_s=SD(y_s)>0,

    beta_s = r_s' z(y_s) / (r_s' h_s).

This measures association with the within-species pattern relative to that species' observed response variation. When sigma_s=0, z(y_s) is undefined. Setting beta_s=0 is not the original estimator.

There is no unique continuous extension at zero: for any centered unit-variance host vector h independent of the controls, responses y=+epsilon*h and y=-epsilon*h both approach zero as epsilon approaches zero, yet the respective standardized slopes remain +1 and -1. Thus the undefined beta is not merely a rounding error that can be repaired by an arbitrary zero.

### B. Response-scale coefficient, a different exploratory estimand

Keep the same edges and standardized predictors but do not divide y_s by its standard deviation:

    gamma_s = r_s' y_s / (r_s' h_s).

For every species with sigma_s>0, the exact OLS scale identity is

    gamma_s = sigma_s * beta_s.

For a constant response, including the four all-zero responses, r_s' y_s=0 because the control matrix includes an intercept. Hence gamma_s=0 is the uniquely defined sample OLS coefficient under this DIFFERENT estimand. It is not an imputed original beta and is not proof that the true biological mechanism is absent.

The interpretation is change in post-IBD p-distance per one within-species SD increase in historical host-climate mismatch, holding the original self-history and current host-climate controls fixed. Multiplying gamma by 100 expresses the response change in percentage points of p-distance. It is still an association, not a causal effect.

A full 259-species average on this scale would be

    G = (1/259) * sum_s gamma_s.

Each species contributes once, but this target differs from mean standardized beta: species with more observed response variation can make a larger absolute contribution. Between-species differences in COI evolution, demographic history, sampling and absolute diversity remain interpretation issues. Changing to G cannot retroactively rescue the original confirmatory result.

## Useful next descriptive analysis, without pretending it is confirmation

1. Keep four distinct response states in the ledger: finite standardized response; exactly invariant response; nonfinite/invalid response; missing input. Never equate missingness with zero variation. For the recorded result, the first two counts are 255 and four.
2. Summarize the already-computed 255 finite original betas descriptively: mean, median, quantiles and counts above/below zero, all explicitly conditional on nonzero observed response variation. Do not present their summary as the original 259-species primary test.
3. Separately summarize the four observed invariant species. Their absence of sampled COI differences is limited to these records and comparable sites; it does not establish genome-wide homogeneity or unrestricted migration.
4. Only after reading the saved coefficient table's actual schema, determine whether it includes each species' unstandardized post-IBD response SD. If both beta_s and sigma_s are present, gamma_s can be recovered by the exact scale identity; otherwise beta_s alone does not contain enough information. Do not invent SDs or infer them from the previous power result.
5. An exploratory G or conditional 255-species p-value would require a distinct, honestly post-outcome analysis and calibration matching its response scale, species set and invariance handling. The old 259-standardized-response reference distribution and 97% simulated power cannot be reused as validation of either new target. No such new empirical test was performed for this note.

This separates the biological questions: how much sampled genetic differentiation exists, and how the variation that exists aligns with historical host-resource opportunity. With only four invariant observations, no fitted model of the probability of invariance is justified by this note; the first component is descriptive, not an inferred two-part ecological mechanism.

## What actually ran in the current continuation

The local Python/container tools returned ClientError and the ZIP content could not be inspected through their runtime. The Files text reader does not expose the CSVs inside this ZIP. Therefore the 255-species coefficient mean, median and sign count have NOT been recomputed here; no value is reported for them, and no replacement primary p-value was generated.

The previously failed persistent save was repaired using the exact exported file snapshots, without depending on the expired container session. The private Library now contains:

- `/TTF/historical_host_memory/one_shot_259_20261012/TTF_hhm259_one_shot_audit_20261012.zip` (59,310 bytes)
- `/TTF/historical_host_memory/one_shot_259_20261012/empirical_result.json`

These were copied by the Files service and their presence was checked by Library listing. The archive's recorded earlier checksum is `d236cb0870f29b7dd3038b2d1b8e338da80a223ed1f57b309f8a0e86a0172123`; a new independent runtime recomputation of that archive checksum was not possible in this turn. Original nucleotide files were not opened or published.

The separate test file `tests/test_hhm_postopening_response_scale_math.py` checks the mathematical claims using synthetic arrays only: zero response-scale slopes for constant response; gamma=SD(y)*beta; partial projection versus full OLS; opposite standardized limits at zero; differing aggregate estimands; and unchanged hash/status of the frozen empirical result. Passing those tests is mathematical/software verification, not a numerical result for the 255 real coefficients and not statistical power validation.

## Sources

- `benchmarks/frozen/historical_host_memory_259_one_shot_empirical_result_v0.1.json`
- `benchmarks/frozen/historical_host_memory_259_zero_response_diagnostic_v0.1.json`
- `docs/supporting/historical_host_memory_qualification_rule_v0.2.json`
- `docs/HISTORICAL_HOST_MEMORY_259_EMPIRICAL_OPENING_V01.md`

All new work in this note is labelled post-outcome. The original frozen record and development-panel closure remain unchanged.
