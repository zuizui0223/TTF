# Lepidoptera trait-gradient study — result spine v0.1

## Central question

Within Lepidoptera, does ecological similarity make post-IBD spatial genetic differentiation more reusable across species once geographic sampling geometry is controlled?

## Design ladder

1. The completed all-animal study provides the broad-domain reference: cross-species post-IBD transfer was not detected in the fresh phylogatR panel.
2. The new study restricts the domain to Lepidoptera and excludes every species from the prior identity-opened Phase-1 route.
3. Within this species-disjoint Lepidoptera domain, the primary test asks whether pairwise transfer congruence increases with a frozen four-axis ecological similarity score (wing size, voltinism, host breadth and habitat affinity), after response-blind geometry control.

The initial estimator v0.1 was retired before empirical opening because it leaked under a geometry-confounded synthetic trap. The successor v0.2 used target-fixed-effect geometry residualization and worst-case nuisance-null envelope calibration. On the exact 240-species character-mask survivor geometry, v0.2 passed all three predeclared gates:

- private-null Type-I: 10/500 rejections; Wilson 95% upper = 0.0364;
- geometry-confounded trap Type-I: 18/500; Wilson 95% upper = 0.0562;
- trait-gradient positive power: 435/500; Wilson 95% lower = 0.8377.

Only after this exact-geometry requalification passed was nucleotide identity opened.

## Primary empirical result

The one-shot primary statistic was -0.0166475 with worst-case envelope p = 0.89 across 110 finite evaluation-target correlations. Fifty target correlations were positive and 60 were negative; the median target correlation was -0.0224.

Frozen decision: **NO_DETECTED_POSITIVE_TRAIT_SIMILARITY_GRADIENT**.

Interpretation: within the exact frozen Lepidoptera survivor domain, ecologically more similar species did not show greater post-IBD spatial-genetic congruence after the frozen geometry control. This does not establish that the true ecological gradient is exactly zero or that ecological traits are universally irrelevant.

## Descriptive second rung

The predeclared within-Lepidoptera geography-conditioned baseline was recovered after the primary decision using the same frozen post-IBD edge responses. Its mean held-out Spearman transfer score was 0.0516743 across 110 supported evaluation species (median species score 0.0551; 62 positive, 48 negative).

This baseline is **descriptive only**. It must not be tested against, or interpreted as significantly larger than, the all-animal result because the taxonomic domain, species panel and survivor geometry differ.

## Exploratory: is transferability itself a species property?

After the frozen primary decision was complete, we asked a separate exploratory question: perhaps ecological similarity fails because transferability is not a property of pairwise similarity but a repeatable property of individual species. We therefore decomposed the exact primary pairwise congruence surface into a source-species **exportability** effect, a target-species **receptivity** effect, and pair-specific residual variation after the same response-blind geometry covariates were removed.

This analysis used 2,428 supported target-source pairs, 106 training species that entered at least one supported pair, and the 110 supported evaluation targets. Source identity accounted for 0.48% of geometry-adjusted pairwise variation, target identity for 1.36%, and 98.16% remained pair-specific/residual. In deterministic 10-fold held-out-pair prediction, adding source and target species identity did not improve RMSE over geometry alone (RMSE 0.265592 versus 0.265568), although the prediction correlation increased only slightly (0.0669 versus 0.0582).

Because this estimand was defined after the primary empirical result, it carries no confirmatory p-value and cannot alter the frozen primary conclusion. Its role is hypothesis-generating. The pattern suggests that transferability is not strongly repeatable as a coarse species-level attribute in this panel. Together with the failed ecological-similarity gradient, this points toward a more relational view: whether spatial genetic structure transfers may depend chiefly on the specific source-target pairing, shared historical exposure, barrier context, or finer interactions rather than on one species' general tendency to export or receive a spatial pattern.

## Ecological message

The useful result is not merely “traits were nonsignificant.” The method was explicitly challenged on the exact deployment geometry and could recover a predeclared trait-gradient world while rejecting geometry-only pseudo-gradients. Yet the empirical Lepidoptera panel showed no positive gradient.

The bounded ecological conclusion is therefore:

> Broad taxonomic restriction can define a more homogeneous comparison domain, but similarity in wing size, voltinism, host breadth and habitat affinity does not predict which Lepidoptera share reusable spatial genetic structure after geography is controlled.

This shifts the biological explanation away from simple contemporary ecological similarity and toward unmeasured lineage-specific history, barriers, demographic responses, or finer biological mechanisms. Those alternatives are hypotheses for independent future data, not post hoc explanations tested here.

## Claim boundaries

- No post-result retuning of traits, geometry covariates, alignment rules, thresholds or nuisance families.
- No claim that the ecological gradient is exactly zero.
- No inferential comparison between the descriptive within-Lepidoptera baseline and the earlier all-animal statistic.
- No claim that any specific historical mechanism caused the observed idiosyncrasy.
- The empirical conclusion is restricted to the frozen 240-species survivor domain and the COI-family response definition.
