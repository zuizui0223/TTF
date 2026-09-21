# Butterfly trait-gradient study with independent Lepidoptera host-resource follow-up — result spine v0.1

## Central question

Within butterflies, does ecological similarity make post-IBD spatial genetic differentiation more reusable across species once geographic sampling geometry is controlled, and does an independent broader Lepidoptera panel reveal a role for native larval host-resource geography?

## Design ladder

1. The completed all-animal study provides the broad-domain reference: cross-species post-IBD transfer was not detected in the fresh phylogatR panel.
2. The trait study restricts the domain to the butterfly species represented in the frozen LepTraits complete-case panel and excludes every species from the prior identity-opened Phase-1 route.
3. Within this species-disjoint butterfly domain, the primary test asks whether pairwise transfer congruence increases with a frozen four-axis ecological similarity score (wing size, voltinism, host breadth and habitat affinity), after response-blind geometry control.

The initial estimator v0.1 was retired before empirical opening because it leaked under a geometry-confounded synthetic trap. The successor v0.2 used target-fixed-effect geometry residualization and worst-case nuisance-null envelope calibration. On the exact 240-species character-mask survivor geometry, v0.2 passed all three predeclared gates:

- private-null Type-I: 10/500 rejections; Wilson 95% upper = 0.0364;
- geometry-confounded trap Type-I: 18/500; Wilson 95% upper = 0.0562;
- trait-gradient positive power: 435/500; Wilson 95% lower = 0.8377.

Only after this exact-geometry requalification passed was nucleotide identity opened.

## Taxonomic, geographic and trait scope of the butterfly panel

A post-result descriptive audit reconstructed the exact frozen panel without changing any inference. The 339 pre-mask species belonged to six butterfly families: Nymphalidae (143), Lycaenidae (64), Hesperiidae (57), Pieridae (56), Papilionidae (18), and Riodinidae (1). After the character-mask gate, the 240 survivors retained the same butterfly-only scope. The 110 geographically supported evaluation species comprised Nymphalidae (59), Hesperiidae (18), Lycaenidae (15), Pieridae (13), and Papilionidae (5). Thus the completed trait-gradient result should not be described as representative of Lepidoptera as a whole.

The evaluation panel was also geographically concentrated. Of 5,456 exact evaluation localities, 2,923 fell in Europe and 1,875 in North America under Natural Earth 110-m country polygons; 357 small-island or ocean localities were unassigned at that map resolution. By dominant continent, 67 evaluation species were North American and 43 European. The median locality latitude was 44.0 degrees (IQR 38.8–48.5 degrees), so the empirical trait-gradient domain is primarily a northern-temperate butterfly comparison.

The null result was not caused by an absence of measured trait variation. Among the 110 supported evaluation species, the frozen wing-size proxy ranged from 1.95 to 11.32, host breadth from 1 to 11 host-plant families, and voltinism included all three frozen categories (U = 28, B = 16, M = 66). Wing size and log host breadth were only weakly correlated (r = 0.099). Using median splits, all four wing-size × host-breadth combinations were represented. These are descriptive scope checks only and do not alter the primary statistic.

Character-mask dropout was reported rather than corrected post hoc. Family-specific dropout fractions were 0.245 for Nymphalidae, 0.281 for Lycaenidae, 0.333 for Hesperiidae, 0.375 for Pieridae, and 0.333 for Papilionidae. No family-specific reweighting or species replacement was introduced after observing these differences.

## Primary empirical result

The one-shot primary statistic was -0.0166475 with worst-case envelope p = 0.89 across 110 finite evaluation-target correlations. Fifty target correlations were positive and 60 were negative; the median target correlation was -0.0224.

Frozen decision: **NO_DETECTED_POSITIVE_TRAIT_SIMILARITY_GRADIENT**.

Interpretation: within the exact frozen butterfly survivor domain, ecologically more similar species did not show greater post-IBD spatial-genetic congruence after the frozen geometry control. This does not establish that the true ecological gradient is exactly zero or that ecological traits are universally irrelevant.

## Descriptive second rung

The predeclared within-butterfly geography-conditioned baseline was recovered after the primary decision using the same frozen post-IBD edge responses. Its mean held-out Spearman transfer score was 0.0516743 across 110 supported evaluation species (median species score 0.0551; 62 positive, 48 negative).

This baseline is **descriptive only**. It must not be tested against, or interpreted as significantly larger than, the all-animal result because the taxonomic domain, species panel and survivor geometry differ.

## Exploratory: is transferability itself a species property?

After the frozen primary decision was complete, we asked a separate exploratory question: perhaps ecological similarity fails because transferability is not a property of pairwise similarity but a repeatable property of individual species. We therefore decomposed the exact primary pairwise congruence surface into a source-species **exportability** effect, a target-species **receptivity** effect, and pair-specific residual variation after the same response-blind geometry covariates were removed.

This analysis used 2,428 supported target-source pairs, 106 training species that entered at least one supported pair, and the 110 supported evaluation targets. Source identity accounted for 0.48% of geometry-adjusted pairwise variation, target identity for 1.36%, and 98.16% remained pair-specific/residual. In deterministic 10-fold held-out-pair prediction, adding source and target species identity did not improve RMSE over geometry alone (RMSE 0.265592 versus 0.265568), although the prediction correlation increased only slightly (0.0669 versus 0.0582).

Because this estimand was defined after the primary empirical result, it carries no confirmatory p-value and cannot alter the frozen primary conclusion. Its role is hypothesis-generating. The pattern suggests that transferability is not strongly repeatable as a coarse species-level attribute in this panel. Together with the failed ecological-similarity gradient, this points toward a more relational view: whether spatial genetic structure transfers may depend chiefly on the specific source-target pairing, shared historical exposure, barrier context, or finer interactions rather than on one species' general tendency to export or receive a spatial pattern.

## Exploratory: which pairwise ecological axis carries the remaining relational signal?

Because the primary composite similarity score was null, we next asked whether that result could reflect cancellation among individual ecological axes. This was defined only after the frozen primary decision and is therefore exploratory. Using the same 2,428 supported target-source pairs, we compared deterministic 10-fold held-out-pair prediction from the frozen geometry covariates alone against geometry augmented by wing-size similarity, host-breadth similarity, voltinism similarity, habitat similarity, the four axes entered separately, and a wing-by-host joint-match term.

None improved held-out-pair RMSE. The geometry-only RMSE was 0.265522. Augmented RMSEs were 0.265603 for wing-size similarity, 0.265701 for host-breadth similarity, 0.265566 for voltinism similarity, 0.265706 for habitat similarity, 0.265782 for wing plus host breadth, 0.265629 for their joint-match product, and 0.266007 when all four axes were entered separately. Prediction correlations likewise did not improve.

This rules out a simple masking explanation within the measured LepTraits variables: the composite null was not produced by one positive trait axis being cancelled by another negative axis. Importantly, host breadth is not the same quantity as the geography of larval host resources. A future independent test can therefore target host-resource geography directly without reinterpreting the present result.

## Independent prospective test: native larval host-resource geography

The completed trait study left one biologically specific alternative unresolved. Host breadth counts how many host taxa a larva can use, but it does not describe **where those resources occur**. We therefore defined a new test on a fully species-disjoint Lepidoptera panel before opening any genetic outcome from that panel.

Larval hosts were taken from the frozen HOSTS database and normalized exactly against WCVP version 13. For each insect species, the primary host-resource footprint was the union of native, extant, non-doubtful WGSRPD level-3 units occupied by its WCVP-resolved larval host plants. Pairwise resource similarity was the Jaccard overlap of these frozen geographic footprints. The new panel excluded all 250 species from the earlier broad-animal identity-opening route, all 339 species from the completed trait-gradient pre-mask universe, and the separately recorded operator-exposure species.

Response-blind screening yielded 500 frozen species before the character-mask gate. The exact survivor panel contained 431 species, inheriting 213 training and 218 evaluation roles; 204 evaluation targets retained the unchanged 500-km geographic support rule, giving 12,358 target-source pairs.

The first host-resource estimator, v0.1, controlled all three nuisance Type-I cells but failed its predeclared positive-power gate: only 115/500 true host-resource-gradient worlds rejected, with Wilson 95% lower power bound 0.195. It was therefore retired before nucleotide identity was opened. Response-blind diagnosis showed that raw host-resource Jaccard still carried strong pairwise resource-breadth similarity. The single v0.2 successor change added pairwise resource-breadth similarity to the frozen nuisance residualization; all other predictor, geometry, alignment, simulator and decision components were unchanged, and a fully disjoint formal seed namespace was used.

On the exact 431-species survivor geometry, v0.2 passed all four gates. Private structure rejected in 1/500 worlds (Wilson 95% upper = 0.0112), the insect-geometry trap in 12/500 (upper = 0.0415), and the host-breadth trap in 29/500 (upper = 0.0821). The true host-resource-gradient cell rejected in 465/500 worlds, with Wilson 95% lower power bound = 0.9042. Only after these four gates passed was nucleotide identity opened once.

The one-shot empirical host-resource statistic was **-0.00704**, with worst-case nuisance-envelope **p = 0.698**. Of 201 finite evaluation-target correlations, 101 were positive and 100 negative; the median was 0.00171 (Q1 = -0.0902, Q3 = 0.0806).

Frozen decision: **NO_DETECTED_POSITIVE_HOST_RESOURCE_GEOGRAPHY_GRADIENT**.

Thus, even when larval resource constraint was represented by the realized native geography of recorded host plants rather than by host breadth alone, species sharing more host-resource geography did not show greater congruence in post-IBD mitochondrial spatial differentiation after insect sampling geometry and simple resource breadth were controlled. This is not evidence that larval resources never affect gene flow; it is evidence that the specific cross-species geographic-overlap gradient tested here was absent in an independently qualified panel.

## Ecological message

The useful result is not merely “traits were nonsignificant.” The method was explicitly challenged on the exact deployment geometry and could recover a predeclared trait-gradient world while rejecting geometry-only pseudo-gradients. Yet the empirical Lepidoptera panel showed no positive gradient.

The bounded ecological conclusion is therefore:

> In the butterfly trait panel, broad ecological similarity did not identify transferable pairs; in a separate, moth-rich broader Lepidoptera panel, overlap in native larval host-resource geography also failed to identify them.

The independent host-resource follow-up strengthens the boundary around the null result: the failure of the coarse LepTraits gradient cannot be attributed simply to host breadth being an inadequate proxy for where larval resources occur. Together with the exploratory finding that source-species and target-species identity explain little pairwise transfer variation, the remaining structure appears predominantly relational or context-specific. Historical co-exposure, barrier-specific responses, demographic history, finer host configuration, symbiont-associated mitochondrial history, and other lineage-pair mechanisms remain hypotheses for independent future data rather than explanations demonstrated here.

## Claim boundaries

- No post-result retuning of traits, geometry covariates, alignment rules, thresholds or nuisance families.
- No claim that the ecological gradient is exactly zero.
- No inferential comparison between the descriptive within-Lepidoptera baseline and the earlier all-animal statistic.
- No claim that any specific historical mechanism caused the observed idiosyncrasy.
- The original trait-gradient conclusion is restricted to its frozen 240-species survivor domain; the independent host-resource conclusion is restricted to its distinct frozen 431-species survivor domain. Both use the declared COI-family response definition.
- The host-resource analysis tests overlap of native host-resource footprints, not a direct causal effect of host-plant distribution on gene flow.
- No post-result tuning of WCVP matching, host footprints, resource-breadth control, geometry covariates or nuisance references is permitted.
