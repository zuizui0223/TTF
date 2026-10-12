# Historical host-memory 259: empirical opening and terminal undefined response

## Current state (2026-10-12)

**The confirmatory nucleotide identities have now been opened.** Older statements that this 259-species panel is unopened are historical, not current.

The user requested continuation after the one-shot review. Authorization was recorded in `benchmarks/frozen/historical_host_memory_259_one_shot_authorization_v0.1.json`, commit `a530d7f58fb798743a776a8b49436ad38d6bfa43`, before accessing the empirical identities. The frozen 259 species and 20,362 original graph edges were processed in the local/private runtime. The original 320-species development identity remains closed. No original sequences, original ZIP, or per-edge genetic distances were published or serialized by the empirical runner.

**Result: `NOT_EVALUABLE_EMPIRICAL_RESPONSE`. The primary statistic and Monte Carlo p-value are undefined and were not computed. This is NOT a biological null and NOT a positive host-memory result.**

## Why it stopped

The primary model requires a standardized response, z(post-IBD turnover), for every species in the fixed 259-species equal-weight panel. Four species have zero uncorrected p-distance on every frozen sampled edge, and hence zero post-IBD response standard deviation:

| Species | Mapped sequences | Localities | Frozen edges | Variable canonical columns |
| --- | ---: | ---: | ---: | ---: |
| Autographa bractea | 16 | 13 | 18 | 0 |
| Bungalotis erythus | 34 | 26 | 64 | 0 |
| Chloantha hyperici | 16 | 14 | 20 | 0 |
| Pennisetia hylaeiformis | 13 | 12 | 15 | 0 |

Their preopened validity-mask hashes match. All original edges have valid comparisons. A separate element-by-element check found zero mismatches across all 185 qualifying sequence-pair comparisons on those 117 edges. This is not a missing-file problem, a failed alignment-overlap gate, or a rounding problem.

With standard deviation zero, z(Y)=(Y-mean(Y))/sd(Y) is undefined. Assigning beta_host=0 by convention or deleting these four species would change the frozen estimand. Neither was done. No replacement 255-species aggregate test or p-value was calculated. The remaining 255 finite species coefficients were produced as part of the authorized one-shot computation and retained locally, but not used to rescue the primary test.

The finding is limited to the sampled COI records and comparable positions. It does not prove genome-wide homogeneity, unrestricted migration, or absence of a historical host effect in those species.

## Why earlier qualification did not detect this

Canonical-site masks certify comparability of positions, not segregating genetic variation. The external predictor-information gate concerns X, not empirical Y. The synthetic world generates continuous private spatial fields and nonzero independent edge noise, so its simulated responses do not exercise an exactly invariant empirical COI response. Its earlier 97% power claim is conditional on that configured simulation family; it did not establish that all actual standardized responses would be defined.

Thus the earlier claim that every prerequisite for the empirical test was satisfied was too strong. Source availability, masks, predictor information and the declared synthetic family passed; empirical response definability had not yet been established.

## Reproducibility and unchanged decision boundary

Before identity opening, the saved private simulator `independent_survivor_synthetic_replay.py` (SHA256 `9c7c2cbe3d56f442215e72887a62ba5019185aff8546c1792f2bfe3230872cc9`) was replayed without parameter changes. It reproduced 30/500 null rejections and 485/500 positive-world rejections. Its 1,999 null-reference statistics were saved before the empirical opening, SHA256 `cbc07fc3d16a601f536acb0b80a922510f82a45f9e2c5f2d243959346f596dab`.

That saved private script uses the canonical 64-hex source ID. The public prototype retained an earlier truncated ID in its seed key; the authorization explicitly binds the actually assessed private replay and its reference series rather than silently treating different seeds as identical.

The one-shot code hashes were frozen in the authorization. Thirteen tests passed before opening. Subsequent result/root-cause checks bring the local suite to sixteen passing tests. All 255 finite coefficients agreed with direct full OLS to within 4.45e-14, and 765 predetermined IBD residual checks agreed with explicit endpoint-disjoint least squares to within 1.12e-16. These are numerical verification checks, not additional ecological hypothesis tests.

Primary result: `benchmarks/frozen/historical_host_memory_259_one_shot_empirical_result_v0.1.json` (SHA256 `5d994a76f7f20a778f0394f545c08a729182ce32298503c128d95c00da29f965`). The GitHub blob matches the local terminal result byte-for-byte.

Root-cause summary: `benchmarks/frozen/historical_host_memory_259_zero_response_diagnostic_v0.1.json`. This is a descriptive replay of the already opened four records and does not alter the primary result.

The execution code, tests, authorization, fixed synthetic reference and terminal records are supplied in the conversation audit bundle; it contains no original nucleotide sequences or raw ZIP. Additional Library uploads failed with `container_session_expired`, so no new persistent Library save is claimed. The original source and earlier saved Library artifacts were not changed.

## Prohibited reinterpretations

Do not call this p=1, p>0.05, evidence for no effect, a failed power test, or a successful 255-species confirmatory test. Do not reset the panel to unopened. Any future redefinition of the response or population must disclose this outcome exposure and cannot retrospectively replace this frozen one-shot result.
