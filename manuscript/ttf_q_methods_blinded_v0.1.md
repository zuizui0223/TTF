# Evaluable envelopes for dyadic ecological predictors: response-blind decomposition of information, calibration and detectability

**TTF-Q methods manuscript v0.1**

Target article type: Research Article, Methods in Ecology and Evolution

## Abstract

1. Pairwise ecological predictors are increasingly used to explain biological relationships between species, populations or sites. Yet a predictor can fail to support inference for different reasons: most of its variation may be endpoint identity, it may be redundant with declared controls, its residual signal may be concentrated in few endpoints, or the intended estimator may be miscalibrated or underpowered. These states are often collapsed into one significance test or power gate.

2. We introduce TTF-Q, a response-blind method that characterizes a dyadic predictor before the biological response is opened. TTF-Q decomposes relation variation into exact source/target fixed-effect survival and nested nuisance-control survival, quantifies residual signal concentration, evaluates null behavior under declared nuisance regimes, and reports minimum detectable effects only where the inferential procedure satisfies a predeclared Type-I qualification ceiling. The joint output is an evaluable envelope rather than a scalar PASS/FAIL label.

3. In pre-result-frozen known-truth simulations, TTF-Q recovered prescribed endpoint, nuisance-control and geographic information-survival fractions across 243 dyadic designs with maximum absolute error of 3.33 × 10^-16. Deliberate source or target concentration was identified in every benchmark comparison. In a second benchmark, four designs failed the same legacy-style beta = 0.03, A = 3 cell but represented distinct states: broad, endpoint-loss and control-redundancy designs had different calibrated detectability envelopes, whereas a source-concentrated design retained the same identifiable relation information as the broad reference but violated the null-qualification ceiling and had no evaluable MDE.

4. We applied TTF-Q to response-blind present and late-Quaternary climate relations across 653 animal species. Historical climate-displacement similarity retained 36–38% of total variation after endpoint identity, present climate, lineage and sampling imbalance; 82–83% of that unique information remained after geographic co-opportunity. On 15,625 development dyads, both climate relations passed null qualification across the nuisance grid. In a separately frozen Palearctic non-Lepidoptera application, 21 species yielded 110 dyads, but LGM similarity retained only 13.1% unique variation and exceeded the predeclared endpoint-concentration ceiling. TTF-Q closed that route before subgroup nucleotide identity or genetic transfer scores were opened.

## Data and code for peer review

An anonymized peer-review archive containing the implementation, machine-readable method contracts, known-truth benchmark rules, frozen benchmark receipts, response-blind climate and prospective-stop receipts, tests and figure-generation code is supplied with the submission. The archive contains no version-control history or author-identifying repository metadata. The manuscript uses no genetic-response quantity from the permanently closed relational v0.3 B/C program.

## Keywords

dyadic data; ecological similarity; evaluable envelope; information decomposition; niche overlap; power analysis; response-blind design; simulation

## 1. Introduction

Ecological inference increasingly relies on pairwise or relational predictors. Examples include environmental similarity between species, geographic co-occurrence, phylogenetic or trait distance, interaction opportunity, landscape connectivity and measures of historical exposure. Such quantities are often entered into models as if the main question were simply whether their fitted coefficient differs from zero.

That framing misses a prior design problem. A source-target relation can appear variable while containing little information that distinguishes one dyad from another after source and target identity are accounted for. It can be almost fully reconstructed from nuisance variables already present in the intended model. Its remaining variation can be concentrated in a small number of endpoints, so a large nominal dyad count gives a misleading picture of where the identifying signal lies. Finally, even a relation with substantial distinct variation can be paired with an inferential procedure that is anti-conservative or weakly powered under the dependence structure of the realized dyadic geometry.

These failure modes are scientifically different. Endpoint redundancy says that apparent pairwise structure is mainly a property of the individual sources or targets. Nuisance redundancy says that the focal relation adds little beyond already-declared predictors. Signal concentration describes where the remaining relation information is located, not whether standard errors are valid. Null miscalibration is a property of the intended inferential procedure under a declared nuisance regime. Low power is meaningful only after calibration has been qualified. Collapsing these states into one PASS/FAIL gate discards information about why a design can or cannot support an intended claim.

Several established methods address pieces of this problem. Ecological modelling has long treated predictor collinearity as a distinct threat to estimation and transferability (Dormann et al. 2013). Ecological niche overlap can be estimated from occurrence data in environmental space, including by kernel-density approaches and Schoener-type overlap measures (Broennimann et al. 2012). Monte Carlo simulation is widely used to estimate power under complex designs (Green & MacLeod 2016), including study-design procedures that map detectable effect sizes from pilot ecological data (Maslen et al. 2023). More generally, simulation has been advocated as a way to understand and validate complex ecological models before empirical use (DiRenzo et al. 2023), and transparent reporting standards for statistical simulation studies now emphasize explicit data-generating mechanisms, complete scenario reporting, Monte Carlo uncertainty and open code (Williams et al. 2024). Multiway and dyadic cluster-robust variance estimators address dependence created when observations share endpoints (Cameron et al. 2011; Aronow et al. 2015). Preregistration and related response-blind practices constrain outcome-dependent analytical choice and clarify when analysis plans change (Gould et al. 2026). None of those components is claimed as new here.

The methodological gap addressed by TTF-Q is how to join them into a single pre-response characterization of a **dyadic predictor itself**. We define an evaluable envelope that keeps five questions distinct:

1. How much focal-relation variation survives exact source and target identity?
2. How much of that surviving relation is distinct from declared nuisance controls?
3. How broadly is the residual relation signal distributed across dyads and endpoints?
4. Under each declared nuisance regime, does the intended inferential procedure satisfy a predeclared Type-I qualification ceiling?
5. Conditional on that qualification, what positive effect sizes are detectably estimable?

TTF-Q is intentionally response blind: all five questions can be answered before the biological outcome being predicted is opened. This ordering is useful not only for confirmatory analyses. It also provides an interpretable diagnostic for exploratory relational predictors, because it separates poor measurement or redundant geometry from the absence of an empirical biological association.

We develop and demonstrate the method in four steps. First, we derive exact nested information-survival identities for fixed dyadic predictor geometry. Second, we validate those quantities and the separation of signal concentration from inferential qualification in known-truth simulations frozen before results were generated. Third, we apply TTF-Q to present-climate and late-Quaternary climate-displacement relations built from occurrence and CHELSA climate data. Fourth, we use the same response-blind logic prospectively to decide whether a separately frozen Palearctic LGM analysis may open its subgroup biological response. No genetic-transferability outcome is used to qualify either ecological application.

## 2. Materials and Methods

### 2.1 Dyadic predictor geometry

Consider a fixed set of directed dyads indexed by source s and target t. Let r be the centered focal relation vector across the n observed dyads. Let F be the exact source/target fixed-effect design matrix and M_F its residual-maker.

The focal relation remaining after endpoint identity is

\[
r_F = M_F r.
\]

Let C denote a matrix of declared nuisance predictors. These may include environmental similarity, geographic opportunity, lineage indicators or sampling-imbalance variables. After applying the same endpoint residualization to the nuisance matrix,

\[
C_F = M_F C,
\]

the focal relation remaining after both endpoint identity and nuisance controls is

\[
r_C = M_{C_F} r_F,
\]

where M_{C_F} is the orthogonal residual-maker for the column space of C_F.

All TTF-Q information quantities are properties of this predictor geometry. They do not use the biological response.

### 2.2 Exact endpoint and control information survival

TTF-Q defines endpoint-identity survival as

\[
S_F
=
\frac{\lVert r_F\rVert^2}{\lVert r\rVert^2},
\]

and control survival conditional on endpoint identity as

\[
S_C
=
\frac{\lVert r_C\rVert^2}{\lVert r_F\rVert^2}.
\]

The total identifiable relation fraction is

\[
S_{\mathrm{total}}
=
\frac{\lVert r_C\rVert^2}{\lVert r\rVert^2}
=
S_F S_C.
\]

Here and throughout, **information** means the squared norm (or sum of
squares) of the standardized focal relation remaining after the declared
projection. It is a predictor-geometry quantity, not Fisher information,
Shannon information, mutual information, or information about the unopened
biological response.

The equality is algebraic. Its value is diagnostic localization rather than novelty of the identity itself. A low S_F means that the focal relation is largely endpoint identity. A low S_C with a moderate or high S_F means that the relation varies within endpoints but is largely reproduced by the declared nuisance-control space.

For nested control sets C_1 and C_2 = [C_1, G] fit on exactly the same dyads, let r_1 and r_2 be the respective residual focal relations. Incremental survival is

\[
S_{2\mid1}
=
\frac{\lVert r_2\rVert^2}{\lVert r_1\rVert^2},
\]

so that

\[
S_{\mathrm{total},2}
=
S_{\mathrm{total},1}S_{2\mid1}.
\]

This makes a sequential information-loss analysis possible. In the empirical climate example, G is directed geographic co-opportunity, so 1 - S_{2|1} is the fraction of previously unique historical-climate relation information removed by that additional geographic control under the declared linear projection.

### 2.3 Residual signal concentration

Information amount and information concentration are distinct. After the final declared projection, TTF-Q assigns each dyad signal mass

\[
w_{st}=r_{C,st}^2
\]

and normalizes those masses to shares pi_st.

The inverse-Herfindahl signal count is calculated for dyads and for source- and target-aggregated shares. For example,

\[
N_{\mathrm{source}}^{\mathrm{signal}}
=
\frac{1}{\sum_s \pi_s^2},
\quad
\pi_s=\sum_t\pi_{st}.
\]

These are concentration indices only. We do not interpret them as independent sample size, effective sample size, degrees of freedom or a replacement for dependence-aware variance estimation.

### 2.4 Dependence-aware inferential procedure

The case-study inferential procedure uses exact source and target fixed effects and two-way source/target cluster-robust inference. Multiple dyads share sources or targets, so treating dyads as independent would be inappropriate. Multiway cluster-robust inference provides a general approach for overlapping clustering dimensions (Cameron et al. 2011), and dyadic-data variance estimators make the same shared-member problem explicit (Aronow et al. 2015).

TTF-Q does not claim this estimator as a new contribution. Instead, it asks whether the intended estimator is suitable on the specific predictor geometry and nuisance regime being considered.

### 2.5 Null qualification

For nuisance-heterogeneity regime A, the method simulates beta = 0 worlds under the declared source effects, target effects, dyad noise and private random-slope structure. Let the Monte Carlo rejection probability be p_0(A) with Wilson interval [L_0(A), U_0(A)].

An amplitude is qualification-compatible when

\[
U_0(A)\le\tau,
\]

where tau is a predeclared Type-I qualification ceiling.

For the present implementation, the one-sided inferential test uses alpha_test = 0.025 and the inherited design-qualification ceiling is tau = 0.05. The distinction is deliberate. Passing the qualification rule means that the Wilson upper bound remains below the declared maximum tolerated null rejection; it does not claim that the simulated rejection probability equals exactly 0.025.

### 2.6 Conditional detectability

For positive focal effect beta, let L_beta(A) denote the Wilson lower bound on detection probability. Given target power pi* = 0.80, the raw grid minimum detectable effect is

\[
D_{\mathrm{raw}}(A)
=
\min\{\beta>0:L_\beta(A)\ge\pi^\*\}.
\]

TTF-Q exposes an evaluable MDE only when the same nuisance-amplitude regime passes null qualification:

\[
D_{\mathrm{eval}}(A)
=
\begin{cases}
D_{\mathrm{raw}}(A), & U_0(A)\le\tau,\\
\mathrm{undefined}, & U_0(A)>\tau.
\end{cases}
\]

This ordering prevents high positive-effect rejection rates from being interpreted as useful power when the same estimator is excessively liberal under beta = 0.

### 2.7 The evaluable envelope

The TTF-Q output is the joint object

\[
\mathcal E
=
(
S_F,
S_C,
S_{2\mid1},\ldots,
N_{\mathrm{dyad}}^{\mathrm{signal}},
N_{\mathrm{source}}^{\mathrm{signal}},
N_{\mathrm{target}}^{\mathrm{signal}},
\{K(A),D_{\mathrm{eval}}(A)\}_{A\in\mathcal A}
).
\]

Optional response-blind measurement-repeatability quantities can be appended when the ecological relation can be reconstructed across independent occurrence or measurement resamples.

The envelope is deliberately not collapsed into one scalar score. Its axes answer different design questions and need not move monotonically together.

### 2.8 Known-truth information benchmark

Before benchmark results were generated, we froze a complete-bipartite simulation design with 20 sources and 20 targets. Focal relations were constructed from mutually orthogonal components in the exact two-way-fixed-effect residual subspace plus an endpoint-only component that is exactly absorbed by source and target identity.

Three quantities were prescribed:

- f, endpoint-retained fraction;
- b, baseline-control survival after endpoint absorption;
- g, geographic survival after the baseline control.

The true retained fractions were therefore

\[
S_{\mathrm{C1}}=fb
\]

and

\[
S_{\mathrm{C2}}=fbg.
\]

We crossed three values of f (0.25, 0.50, 0.75), three values of b (0.25, 0.50, 0.80), three values of g (0.25, 0.50, 0.80), three signal-concentration modes (broad, source-concentrated, target-concentrated) and three deterministic seeds, giving 243 predeclared scenarios.

Exact-recovery error and the source/target concentration diagnostics were computed for every scenario. In line with recent recommendations for statistical simulation studies in ecology, the benchmark contract records the data-generating mechanism, target quantities, complete factor grid, deterministic seed rule and all predeclared scenarios rather than selecting favourable cells after simulation (Williams et al. 2024).

### 2.9 Known-truth detectability benchmark

A second benchmark used 30 sources and 30 targets and compared four predeclared predictor geometries:

- **reference:** f = 0.75, b = 0.80, g = 0.80, broad residual signal;
- **endpoint loss:** f = 0.25, b = 0.80, g = 0.80, broad signal;
- **control redundancy:** f = 0.75, b = 0.25, g = 0.80, broad signal;
- **source concentration:** f = 0.75, b = 0.80, g = 0.80, source-concentrated residual signal.

An initial frozen benchmark established that raw MDEs alone were inadequate because some beta = 0 cells exceeded the intended Type-I qualification ceiling. That diagnostic result was retained rather than overwritten. A second contract was then frozen before rerunning the same scientific scenarios with 10,000 beta = 0 worlds per nuisance-amplitude cell and 1,000 positive-effect worlds.

The positive-effect grid was beta = 0.02, 0.03, 0.05, 0.08, 0.10, 0.15 and 0.20. Private heterogeneity was A = 0, 1, 2 and 3. Source and target intercept SD were each 0.18, dyad-noise SD was 0.32, and private random-slope SD was 0.12A. The response transform was tanh.

For comparison with the earlier binary qualification approach, beta = 0.03 at A = 3 was retained as a legacy-style reference cell but had no authority to replace the full envelope.

### 2.10 Response-blind climate case study

#### Species and occurrence domain

The climate example used an independently selected animal species domain from the earlier TTF relational program. Candidate selection and the response firewall had been frozen before any Study-C genetic response. The present TTF-Q analysis reuses only the response-blind ecological assets.

Occurrence records came from GBIF under the pre-existing deterministic filtering rules: records from 2010–2026 with geographic coordinates, duplicate removal, 10-km thinning, a maximum of 200 retained records per species and a minimum of 30 retained records per species.

After climate-raster validity filtering, 653 species retained sufficient ecological data. No species was lost from the admissible set by the three climate-invalid occurrence rows removed in the fresh TTF-Q reconstruction.

#### Historical climate-displacement relation

Historical climate data came from CHELSA-TraCE21k V1.0 at approximately 1-km resolution (Karger et al. 2023; Karger 2025). Four bioclimatic variables were used at 21 ka BP and 0 BP:

- mean annual near-surface air temperature;
- annual temperature range;
- annual precipitation;
- precipitation seasonality.

For each retained occurrence coordinate, the historical displacement vector was the 21-ka value minus the 0-BP value for the four variables. Pooled displacement vectors were standardized, projected to a canonical-sign PCA, reduced to the first two whitened axes, and represented on a fixed 50 × 50 grid using isotropic Gaussian kernel densities. Pairwise historical similarity R_hist was Schoener's D between the resulting species densities.

#### Current climate relation

The current-climate nuisance representation used CHELSA V2.1 1981–2010 bioclimatic variables bio1, bio7, bio12 and bio15 (Brun et al. 2022a,b). The same density and Schoener-overlap machinery was used to obtain R_current.

Schoener-type overlap in gridded environmental space is established in ecological niche-comparison methodology (Broennimann et al. 2012); TTF-Q does not treat the overlap metric itself as a methodological innovation.

#### Panels and dyads

Species were assigned to panels by hashes frozen before response opening. The development panel contained 250 species split into 125 sources and 125 targets, producing 15,625 directed dyads. The confirmatory ecological panel contained 403 remaining species and 40,602 directed dyads. Both panels are reported; neither was selected based on the TTF-Q result.

#### Baseline controls and geography

The baseline response-blind controls were same-class, same-order and same-family indicators plus the absolute log source/target retained-occurrence-count ratio. Historical climate was additionally conditioned on current-climate similarity.

Geographic co-opportunity was defined continuously for each directed source-target dyad as the fraction of target occurrence coordinates whose nearest source occurrence was within 500 km great-circle distance. No hard geographic filter was applied in TTF-Q.

### 2.11 Prospective stopping application

To test whether the evaluable-envelope logic could govern a real analysis rather than only characterize completed ecological geometry, we prospectively froze a secondary Palearctic insect application before constructing its subgroup genetic response. The parent 211-species aggregate genetic result was historically known, so this was not treated as an independent confirmatory dataset. However, species-level subgroup genetic scores, pairwise transfer scores, the subgroup focal coefficient and subgroup nucleotide identity were not used to define or qualify the new design.

The formal panel contained 23 terrestrial Palearctic non-Lepidoptera insects assigned response-blindly to 11 source and 12 target species with no role overlap. External occurrence rules, current and 21-ka climate assets, climate variables, LGM relation, controls and TTF-Q thresholds were frozen before the focal relation was evaluated. After the minimum-30 occurrence rule, 21 species remained: 10 sources, 11 targets and 110 directed dyads.

The focal relation, R_LGM, was Schoener's D between normalized continuous 21-ka climatic-opportunity surfaces estimated from the frozen four-variable CHELSA climate space. R_present, directed geographic co-opportunity, lineage indicators, occurrence-count imbalance and exact source/target fixed effects were declared controls. The prospective opening rule required at least 15% total unique focal-relation variation, at least 50% effective source and target signal breadth, and no single source or target carrying more than 15% of residual signal. A response-blind execution amendment, frozen after relation construction but before any genetic response or TTF-Q outcome was opened, allowed deterministic necessary gates to be evaluated first and the analysis to stop immediately if any failed. No threshold was changed.

The application was defined to yield either authorization for one pre-specified genetic analysis or `NOT_EVALUABLE`. This binary authorization was an application-specific governance rule layered on the multidimensional TTF-Q envelope; it did not replace the envelope with a universal score or threshold. A qualification failure could not be interpreted as a biological null and could not be followed by a replacement subgroup, historical predictor or retuned threshold.

### 2.12 Response firewall and analytic provenance

TTF-Q was developed after the relational v0.3 B/C program had closed with no evaluable genetic test. That historical terminal state is retained unchanged. The fresh method program did not inspect nucleotide identity, pairwise genetic distances, T_st or the genetic relational coefficients.

Machine-readable contracts were committed before each known-truth benchmark stage. When the first detectability benchmark revealed that raw MDEs could be reported despite unacceptable null behavior, the result was preserved as a diagnostic audit. A revised calibration rule and higher-precision null simulation were frozen before the second benchmark was run. This follows the broader principle that planned and post-result analytical changes should remain distinguishable (Gould et al. 2026).

### 2.13 AI-assisted development disclosure

OpenAI ChatGPT (GPT-5.6 Sol) was used during method development and manuscript preparation to assist with code drafting, test generation, documentation, repository auditing and language editing. AI-generated suggestions were not treated as scientific evidence. Code supporting reported results was executed through the repository test and workflow suite, and reported scalar results were locked to machine-readable frozen receipts. The authors remain responsible for the design, code, analyses, interpretation and manuscript content.

## 3. Results

### 3.1 TTF-Q exactly recovered known information losses

Across all 243 known-truth decomposition scenarios, the maximum absolute error across endpoint survival, baseline-control survival, geographic survival, C1 total unique fraction and C2 total unique fraction was

\[
3.33\times10^{-16}.
\]

Thus the implementation recovered the prescribed nested information geometry to floating-point precision. This first benchmark is primarily an implementation verification of the exact projection algebra; it is not, by itself, evidence that TTF-Q outperforms alternative design diagnostics. The non-trivial inferential validation comes from the second benchmark, where information amount, concentration, null qualification and detectability are deliberately allowed to diverge.

Signal concentration was also identified independently of total information. Source-concentrated designs had lower effective-source signal counts than their broad counterparts in 100% of predeclared comparisons; target-concentrated designs had lower effective-target counts in 100%. Median effective endpoint counts under deliberate concentration were approximately 22–23% of those in matched broad-signal designs (Fig. 4).

### 3.2 One binary gate conflated distinct known-truth states

All four known-truth detectability designs received the same negative result at the legacy-style beta = 0.03, A = 3 reference cell.

The calibrated evaluable envelopes were nevertheless different (Fig. 2; Table 1).

**Table 1. Known-truth calibration and evaluable minimum-detectable-effect benchmark.**

| Constructed design | Total unique relation fraction | Null qualification | Evaluable MDE at A = 0 / 1 / 2 / 3 |
| --- | ---: | --- | --- |
| Reference | 0.48 | pass at all A | 0.05 / 0.08 / 0.08 / 0.10 |
| Endpoint loss | 0.16 | pass at all A | 0.08 / 0.10 / 0.15 / 0.20 |
| Control redundancy | 0.15 | pass at all A | 0.08 / 0.10 / 0.15 / 0.20 |
| Source concentration | 0.48 | fail at all A | undefined at all A |

The reference, endpoint-loss and control-redundancy designs therefore differed in detectability despite sharing the same binary label. Endpoint loss and control redundancy had the same MDE vector but were distinguished by where their information was removed: endpoint identity versus nuisance-control redundancy.

The source-concentrated design was qualitatively different. It retained the same total unique relation fraction as the reference design, but its effective-source signal count fell from 28.4 to 6.15. Its beta = 0 rejection rate was 6.1–6.8% across nuisance amplitudes, with Wilson upper bounds above tau = 0.05 in every case. Raw power values were therefore not licensed as evaluable MDEs.

Among the six pairwise comparisons of the four designs, five had the same binary label but different evaluable-MDE vectors. The remaining endpoint-loss versus control-redundancy pair had the same MDE vector but different information-loss decomposition.

### 3.3 Present and historical climate carried distinct pairwise information

On the response-blind climate data, historical and current relation values were related but far from identical. Spearman correlation between R_hist and R_current was 0.621 in the development panel and 0.592 in the confirmatory ecological panel.

For the current-climate relation, 41.6% of standardized relation variation survived exact source and target identity in the development panel (Fig. 1). After lineage and sampling-imbalance controls, 40.6% of total variation remained unique. On the shared development geometry after adding continuous geographic co-opportunity, the total unique fraction was 0.285.

Historical climate showed a different pattern. After source and target identity, 67.3% of historical relation variation remained in development and 66.2% in confirmatory. Conditioning additionally on present climate, lineage and sampling imbalance left total unique fractions of 0.365 and 0.381, respectively.

### 3.4 Most history-specific information survived geographic adjustment

Adding continuous 500-km geographic co-opportunity reduced the historical unique fraction from 0.365 to 0.299 in development and from 0.381 to 0.316 in confirmatory.

The corresponding nested geography-survival fractions were 0.820 and 0.830. Thus approximately 17–18% of previously unique historical relation information was attributable to the declared geographic co-opportunity control, while approximately 82–83% remained after that adjustment.

Residual historical signal was not dominated by one or two endpoints. After geographic adjustment, the development panel had inverse-Herfindahl signal counts of 84.8 sources and 78.8 targets, with maximum source and target signal shares of 4.0% and 4.7%. The confirmatory panel had effective counts of 143.5 sources and 145.6 targets, with maximum shares of 2.2% and 2.5%.

### 3.5 Climate-case detectability depended on nuisance heterogeneity

On the shared 15,625-dyad development geometry, all tested nuisance-amplitude levels satisfied the inherited tau = 0.05 null-qualification ceiling for both climate relations.

For current climate, beta = 0 rejection rates across A = 0, 1, 2 and 3 were 0.032, 0.027, 0.032 and 0.018; the maximum Wilson 95% upper bound was 0.0448. The evaluable grid MDE vector was 0.02, 0.02, 0.05 and 0.05 (Fig. 3).

For historical climate given current climate, beta = 0 rejection rates were 0.022, 0.031, 0.023 and 0.033; the maximum Wilson upper bound was 0.0460. The evaluable grid MDE vector was 0.02, 0.03, 0.05 and 0.08.

The old beta = 0.03, A = 3 reference point therefore fell below the evaluable MDE for both relations, even though both designs were qualification-compatible and informative over substantial lower-heterogeneity or larger-effect regions.

### 3.6 TTF-Q prospectively stopped a plausible LGM analysis before the biological response

The Palearctic application passed its external-data gate with 21 species, 10 source clusters, 11 target clusters and 110 directed dyads. Its LGM and present-climate relations were strongly correlated (Pearson r = 0.861). After exact endpoint identity and all declared controls, the total unique fraction of R_LGM was 0.131, below the prospectively frozen minimum of 0.15.

Residual relation information was also too concentrated. The inverse-Herfindahl effective-source and effective-target counts were 7.60 of 10 and 7.68 of 11, respectively, but the largest single-source share was 21.3% and the largest single-target share 19.3%, both above the frozen 15% ceiling. Thus all three deterministic necessary opening gates—unique information, source signal breadth and target signal breadth—failed.

Under the response-blind short-circuit rule, those failures were sufficient to classify the relation as `NOT_EVALUABLE_NONLEPIDOPTERA_PALEARCTIC_LGM_RELATION`. The subgroup biological response remained closed: no subgroup nucleotide identity, pairwise source-to-target transfer scores or focal genetic coefficient were computed. A full synthetic run that had already begun before the short-circuit amendment was frozen later also failed null qualification, but that result is reported only as non-authoritative corroboration of the prospectively determined stop.

## 4. Discussion

TTF-Q replaces a binary question—whether a relational predictor passes one qualification cell—with a structured description of what information the predictor contains and where inference is licensed.

The known-truth simulations demonstrate why this distinction matters. Endpoint loss, nuisance redundancy and concentrated signal can all lead to the same negative binary classification while representing different statistical states. Two of those mechanisms reduced identifiable relation information and worsened calibrated MDEs. Concentration, by contrast, could leave total unique information unchanged while damaging null behavior. No single scalar summary recovered those distinctions.

### 4.1 What is new and what is not

The individual components of TTF-Q are deliberately conventional. Predictor collinearity diagnostics and remedies are well established in ecology (Dormann et al. 2013). Projection on fixed effects and nuisance covariates is linear-model geometry. Multiway and dyadic cluster-robust inference are established approaches (Cameron et al. 2011; Aronow et al. 2015). Simulation-based power and study-design analysis are established (Green & MacLeod 2016; Maslen et al. 2023), and known-truth simulation is already advocated for validating complex ecological models and statistical methods (DiRenzo et al. 2023; Williams et al. 2024). Environmental niche overlap based on occurrence densities is established (Broennimann et al. 2012). Response-blind planning and preregistration principles are also increasingly developed for ecological modelling (Gould et al. 2026).

The contribution is the **joint response-blind method object**. TTF-Q treats predictor information survival, information-loss location, residual endpoint concentration, null qualification and conditional detectability as separate axes of one evaluable envelope. The exact nested identities make the information-loss axes auditable, while known-truth benchmarks show that the axes are not interchangeable.

This distinction is especially important for methods evaluation. A procedure may have low power because the focal relation contains little unique information, because the unique information is obscured by nuisance heterogeneity, or because the inferential procedure is not qualification-compatible under the realized dependence structure. Those states imply different next actions. More ecological observations may improve relation measurement. A different sampling design may reduce endpoint concentration. A different estimand may be needed if nuisance controls absorb nearly all focal information. A different inferential procedure may be needed when null behavior violates the declared ceiling. Simply reporting “underpowered” hides these choices.

### 4.2 Why the climate result is ecologically useful without a biological response

The climate application also yields a substantive result about ecological predictor geometry. Present-day environmental similarity was not merely endpoint identity or geography. Historical climate-displacement similarity in turn contained substantial additional pairwise information after present climate was included.

The historical result is strongest after the explicit geographic decomposition: roughly 82–83% of the history-specific information identified before geography remained after continuous co-opportunity adjustment in two independently assigned panels. This does not show that historical climate predicts genetic transferability. It shows something narrower and cleaner: species pairs that are similar in present climate and spatial opportunity can still differ systematically in the similarity of their late-Quaternary climate-displacement distributions.

This distinction matters for comparative biogeography. Present ecological resemblance can arise through different historical trajectories. A relational variable intended to represent history should therefore be evaluated for information beyond present state and shared geography before it is used as an explanatory covariate.

### 4.3 The same ecological idea can be evaluable in one design and non-evaluable in another

The prospective stopping application adds a distinction that simulations alone cannot demonstrate. Historical-climate relations were informative and qualification-compatible on the broad 653-species climate geometry, yet a biologically related LGM climatic-opportunity predictor failed prospectively in the narrower Palearctic insect design. TTF-Q therefore does not classify a predictor label—such as “historical climate”—as intrinsically good or bad. It evaluates the information geometry and inferential support created by a particular predictor, control set, endpoint composition and dyad structure.

The Palearctic failure was also diagnostically specific. LGM similarity was strongly redundant with present-climate similarity, leaving only 13.1% unique variation after the declared model structure, and the remaining signal was disproportionately carried by a small number of source and target species. Those are design limitations that were knowable before the subgroup genetic response was inspected. The correct action was therefore to stop, not to interpret an eventual weak coefficient, change the historical predictor or search for a favorable subgroup.

This prospective stop illustrates the practical value of response-blind qualification. It can prevent a biologically plausible analysis from becoming a post-outcome search for an estimable effect. The stop rule is downstream governance of outcome opening, not a replacement for the evaluable envelope itself. Equally, `NOT_EVALUABLE` is not evidence against the biological hypothesis: the method declined to ask the outcome question under that realized design.

### 4.4 Relation information is not outcome importance

A high unique fraction does not imply a large biological effect. TTF-Q measures the amount of focal predictor variation available to identify an effect under a declared design; it does not estimate that effect. Likewise, the 82–83% historical geography-survival result is not a variance partition of any organismal response.

This separation is intentional. Predictor geometry is often knowable before the response is opened, which makes it suitable for prospective design and for transparent decisions about whether a proposed relational analysis is worth pursuing.

### 4.5 Signal concentration is not effective sample size

The known-truth concentration experiment provides a caution against overinterpreting inverse-Herfindahl endpoint counts. Source concentration was detected strongly, but its relationship with raw power was not a simple monotone reduction. In the higher-precision benchmark, its main failure was inflated null rejection under the intended estimator, not a uniquely ordered loss of MDE.

We therefore use effective dyad/source/target counts only as labels for signal concentration. Standard errors and qualification remain the responsibility of the dependence-aware inferential procedure and its null simulations.

### 4.6 Qualification is conditional on the declared synthetic worlds

The evaluable envelope is not universal. Null qualification and MDE depend on the simulated nuisance family. Unrepresented forms of private structure can still invalidate inference. The appropriate use of TTF-Q is therefore to make the synthetic world family explicit, test scientifically plausible nuisance regimes and report the envelope conditional on those regimes.

In the present case, the nuisance family was inherited from a prospectively frozen relational design: source and target intercept heterogeneity, dyad noise, and source- and target-specific random slopes whose SD increased with amplitude A. Other applications should define their own nuisance family before outcome inspection.

### 4.7 Measurement repeatability remains an extension

The current paper validates predictor non-redundancy, nested information loss, signal concentration, null qualification and detectability. TTF-Q also implements a response-blind relation-repeatability layer based on repeated reconstruction of fixed dyads, but the empirical climate example did not require full occurrence-bootstrap reconstruction to establish the core method.

Future applications can add measurement uncertainty as an outer layer: occurrence resampling changes the estimated relation, while TTF-Q then evaluates how that propagated measurement uncertainty changes information survival and the evaluable envelope. This is a natural extension rather than a prerequisite for the validated core method.

## 5. Conclusions

A dyadic ecological predictor can fail before any biological response is analysed, but “fail” is not one statistical state.

TTF-Q decomposes that pre-response problem into identifiable relation information, the location of information loss, residual signal concentration, null qualification and conditional detectability. Known-truth simulations show that these axes can be recovered exactly and can distinguish designs that a single binary power gate treats as identical.

In the response-blind 653-species climate case, present and historical climate each carried substantial distinct pairwise information, and most history-specific information survived explicit geographic adjustment. In contrast, a separately frozen 21-species Palearctic LGM application failed the necessary information and endpoint-support gates, and TTF-Q stopped the analysis before its subgroup genetic response was opened.

The practical output of TTF-Q is therefore not permission to declare a predictor “good” or “bad”. It is a map of the inferential domain that exists for a specific realized design before the outcome is opened—and, when that domain is inadequate, a principled reason not to run the biological test.

## Acknowledgements

Omitted for double-anonymous review.

## Author contributions

Omitted for double-anonymous review.

## Conflict of interest

To be completed.

## Data availability

An anonymized peer-review archive containing all non-standard analysis and simulation code, frozen method contracts, scalar result handoffs and tests is supplied with the submission. Public ecological source datasets are identified in the Methods and reference list. Repository identity and permanent archival DOI will be restored after peer review/acceptance in accordance with journal policy.

## References

Aronow, P. M., Samii, C. & Assenova, V. A. (2015). Cluster-robust variance estimation for dyadic data. *Political Analysis*. https://doi.org/10.1093/pan/mpv018

Broennimann, O., Fitzpatrick, M. C., Pearman, P. B., Petitpierre, B., Pellissier, L., Yoccoz, N. G., Thuiller, W., Fortin, M.-J., Randin, C., Zimmermann, N. E., Graham, C. H. & Guisan, A. (2012). Measuring ecological niche overlap from occurrence and spatial environmental data. *Global Ecology and Biogeography*, 21, 481–497. https://doi.org/10.1111/j.1466-8238.2011.00698.x

Cameron, A. C., Gelbach, J. B. & Miller, D. L. (2011). Robust inference with multiway clustering. *Journal of Business & Economic Statistics*, 29, 238–249. https://doi.org/10.1198/jbes.2010.07136

DiRenzo, G. V., Hanks, E. & Miller, D. A. W. (2023). A practical guide to understanding and validating complex models using data simulations. *Methods in Ecology and Evolution*, 14, 203–217. https://doi.org/10.1111/2041-210X.14030

Dormann, C. F., Elith, J., Bacher, S., Buchmann, C., Carl, G., Carré, G., García Marquéz, J. R., Gruber, B., Lafourcade, B., Leitão, P. J., Münkemüller, T., McClean, C., Osborne, P. E., Reineking, B., Schröder, B., Skidmore, A. K., Zurell, D. & Lautenbach, S. (2013). Collinearity: a review of methods to deal with it and a simulation study evaluating their performance. *Ecography*, 36, 27–46. https://doi.org/10.1111/j.1600-0587.2012.07348.x

Gould, E., Jones, C. S., Yen, J. D. L., Fraser, H. S., Wootton, H. F., Good, M. K., Duncan, D. H., Hauser, C. E., Wintle, B. C. & Rumpff, L. (2026). ‘But I can't preregister my research’: Improving the reproducibility and transparency of ecology and conservation with adaptive preregistration for model-based research. *Methods in Ecology and Evolution*, 17, 1768–1787. https://doi.org/10.1111/2041-210X.70311

Green, P. & MacLeod, C. J. (2016). SIMR: an R package for power analysis of generalized linear mixed models by simulation. *Methods in Ecology and Evolution*, 7, 493–498. https://doi.org/10.1111/2041-210X.12504

Brun, P., Zimmermann, N. E., Hari, C., Pellissier, L. & Karger, D. N. (2022b). Global climate-related predictors at kilometer resolution for the past and future. *Earth System Science Data*, 14, 5573–5603. https://doi.org/10.5194/essd-14-5573-2022

Brun, P., Zimmermann, N. E., Hari, C., Pellissier, L. & Karger, D. N. (2022a). CHELSA-BIOCLIM+ A novel set of global climate-related predictors at kilometre-resolution. EnviDat. https://doi.org/10.16904/envidat.332

Karger, D. N., Nobis, M. P., Normand, S., Graham, C. H. & Zimmermann, N. E. (2023). CHELSA-TraCE21k – high-resolution (1 km) downscaled transient temperature and precipitation data since the Last Glacial Maximum. *Climate of the Past*, 19, 439–456. https://doi.org/10.5194/cp-19-439-2023

Karger, D. N. (2025). CHELSA-TraCE21k-centennial-bioclim and topographic data since the Last Glacial Maximum. EnviDat. https://doi.org/10.16904/envidat.691

Maslen, B., Popovic, G., Lim, M., Marzinelli, E. & Warton, D. (2023). How many sites? Methods to assist design decisions when collecting multivariate data in ecology. *Methods in Ecology and Evolution*, 14, 1564–1573. https://doi.org/10.1111/2041-210X.14094


## Figure legends

**Fig. 1. Response-blind relational information survival.** Fraction of total focal-relation variance retained on the shared development geometry after exact source/target identity, the baseline declared controls, and continuous geographic co-opportunity. Current climate and historical climate are shown separately; the historical baseline-control stage additionally conditions on current-climate similarity.

**Fig. 2. Known-truth designs share a binary label but not an evaluable envelope.** Calibration-qualified grid minimum detectable effects across private-heterogeneity amplitudes A = 0–3 for the broad reference, endpoint-loss and control-redundancy designs. The predeclared source-concentrated design is not assigned an evaluable MDE because its beta = 0 cells exceed the declared Type-I qualification ceiling at every tested A.

**Fig. 3. Empirical response-blind climate relation envelopes.** Calibration-qualified grid minimum detectable effects across private-heterogeneity amplitudes for present-climate similarity and historical climate-displacement similarity conditional on present climate on the 15,625-dyad development geometry. No genetic-response quantity is used.

**Fig. 4. Equal unique relation information can have different endpoint concentration and null qualification.** The known-truth reference and source-concentrated designs retain the same total unique focal-relation fraction, but deliberate concentration lowers the inverse-Herfindahl effective source count and is accompanied by failure of the declared null-qualification ceiling under the intended estimator.

**Fig. 5. A prospective application stopped before the biological response was opened.** Response-blind Palearctic LGM climatic-opportunity similarity was evaluated on 110 directed dyads from 21 species. The panel failed the prospectively frozen unique-information and endpoint-concentration gates, so the planned subgroup genetic response remained unopened. The figure contrasts this realized non-evaluable design with the broader response-blind climate geometry without interpreting the stop as a biological null.


Williams, C., Yang, Y., Lagisz, M., Morrison, K., Ricolfi, L., Warton, D. I. & Nakagawa, S. (2024). Transparent reporting items for simulation studies evaluating statistical methods: Foundations for reproducibility and reliability. *Methods in Ecology and Evolution*, 15, 1926–1939. https://doi.org/10.1111/2041-210X.14415
