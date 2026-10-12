# HHM 259: locality-level influence diagnostic

Updated 2026-10-12. This is a post-outcome exploratory diagnostic, not a new preregistration, a repaired confirmatory test or authorization to open the 320 development species.

## Current evidence and provenance

The original 259-species / 20,362-edge confirmatory test remains `NOT_EVALUABLE_EMPIRICAL_RESPONSE`: four invariant responses made its standardized estimand undefined. The later unscaled-response exploration is explicitly a different estimand.

The public aggregate at commit `6d1a6a3c95f67834bb3b9f8d3d4fd97d9f019dd8`, `benchmarks/exploratory/historical_host_memory_259_mean_influence_v0.1.json`, already contains a descriptive response-SD rank split. Its highest rank group comprises 63 species, with 32 positive and 31 negative slopes, and accounts for 0.9858082133271349 of the all-species SIGNED NET unscaled slope total. This is not a fraction of biological variation explained or gross positive contributions. The group was defined after seeing the outcomes. Its response SD is that of the post-IBD residual, not raw COI divergence. The coefficient/response-SD correlation is 0.03273256991231184, so the previously reported 31.2% algebraic covariance contribution is not a strong cross-species correlation or a causal percentage.

The current chat attachment `influence_result.json` and the GitHub aggregate agree on their coefficient input digest `a102adacf4532eacb35d45360e8b4515d5049622c4c6e36c9bc05a5eb5d0f1f2` and their shared mean/tail summaries, but have different schemas, scope hashes and analysis-code hashes. They are not claimed to be the same complete analysis artifact. In particular, the 63-species grouping is verified here from the GitHub record, not newly recomputed in this session.

## Question

Are large partial host-history slopes supported across multiple localities, or are they highly sensitive to one locality and its incident edges? Large post-IBD residual variation can reflect biological variation, sampling geometry and the correction operation; the aggregate summaries alone do not distinguish these mechanisms.

## Exact diagnostic operation

The new module `scripts/historical_host_memory_locality_sensitivity.py` takes ONE species' already-authorized arrays in memory: original integer edge endpoints, chord distances in kilometres, the three original external predictor columns M_host/M_self/D_host_0BP, and original edge p-distances.

For a future private empirical replay:

1. Verify the untouched original archive, candidate, geometry, survivor role, predictor and prior-mask hashes. Use only the original 259 already-opened survivors. Reconstruct the same original edge p-distances privately, without writing sequence or edge-distance files. Verify every full-sample baseline against the earlier source-bound coefficients before interpreting a deletion.
2. Standardize the three predictors ONCE using the full original species graph; keep those units fixed across all deletions. This avoids confusing a change of predictor units with locality influence.
3. For every original locality, temporarily omit all incident edges. Keep all other original edges; never create or rewire edges.
4. Recompute every retained edge's geographic IBD fit using only retained edges sharing NEITHER endpoint of that edge. Require at least five training edges and nonsingular log-distance variation. Do not reuse residuals corrected using the omitted locality.
5. Refit the unscaled post-IBD response on the three full-sample-scaled predictors plus an intercept. Also report the slope divided by the deletion-specific response SD when defined. The latter retains original X units, so it is not a fully restandardized subset model.
6. Preserve every deletion attempt, including NOT_ESTIMABLE cases from sparse geometry or rank deficiency. A stable-sign statement for all deletions is unavailable when coverage is incomplete. The zero-variation species remain zero-slope under the unscaled model and undefined under response standardization.
7. Report coefficient ranges, sign reversals and computability coverage as descriptive diagnostics only. No p-values, confidence intervals, new subgroup tests, reference-null reuse or post-outcome species/edge eligibility changes.

This is a sensitivity calculation, not a leave-one-locality-out prediction evaluation: predictors remain the original frozen external quantities, and the operation is not a demonstration of out-of-sample ecological transfer.

## Implementation versus empirical completion

A standalone in-memory calculation helper and synthetic-only test workflow have been added. The tests compare the fast disjoint regression with explicit fits, check every deletion against full OLS, verify that an omitted locality cannot affect its deletion result through reused IBD residuals, and cover zero response, loss of estimability and invalid graph inputs.

**No actual locality-deletion results were computed in this session.** Both private Python and container execution returned infrastructure `ClientError`, including after the original survivor CSVs were materialized. The public CI tests are artificial-data checks only; neither the private original archive nor empirical edge-distance vectors have been sent to GitHub. The helper itself is not a source-verifying empirical driver.

## Durable preservation

The prior conversation snapshots were successfully copied, without modifying their bytes or requiring the failed local runtime, into the owner's personal Library folder `/TTF/historical_host_memory/exploratory_followup_20261012/`: the exploratory bundle, influence bundle, both aggregate result JSONs, the influence README and the prior GitHub patch. The private derived-coefficient files remain inside owner-only bundles. No new empirical coefficients, locality robustness claims or biological conclusions have been generated.
