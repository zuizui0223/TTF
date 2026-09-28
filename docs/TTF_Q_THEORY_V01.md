# TTF-Q theory v0.1 — nested information survival and calibrated evaluable envelopes

## Scope

TTF-Q characterizes a dyadic predictor before the biological response is opened.
The quantities below are properties of the focal predictor and its intended
design matrix. They are **not** partitions of biological-response variance.

Throughout TTF-Q, the word **information** is operational shorthand for the
squared Euclidean norm (equivalently, sum of squares) of a standardized focal
relation that remains after declared linear projections. It is not Fisher
information, Shannon information, mutual information, or information about an
unopened biological response.

Let \(r \in \mathbb{R}^n\) be the centered focal source-target relation on a
fixed dyad set. Let \(F\) denote the exact source/target fixed-effect design,
and let \(M_F\) be its residual-maker. Let \(C\) be a matrix of declared
nuisance predictors and define

\[
C_F = M_F C.
\]

The focal relation remaining after endpoint identity is

\[
r_F = M_F r.
\]

The relation remaining after endpoint identity and the declared controls is

\[
r_C = M_{C_F} r_F,
\]

where \(M_{C_F}\) is the orthogonal residual-maker for the column space of
\(C_F\).

TTF-Q reports squared-norm survival fractions because the predictor geometry is
fixed before any biological response is opened.

## Proposition 1 — exact two-stage factorization

Define endpoint-identity survival as

\[
S_F = \frac{\lVert r_F\rVert^2}{\lVert r\rVert^2},
\]

and control survival conditional on endpoint identity as

\[
S_C = \frac{\lVert r_C\rVert^2}{\lVert r_F\rVert^2}.
\]

Then total identifiable relational information is

\[
S_{\mathrm{total}}
=
\frac{\lVert r_C\rVert^2}{\lVert r\rVert^2}
=
S_F S_C.
\]

This is an algebraic identity, not an empirical discovery. Its methodological
use is localization of information loss:

- \(1-S_F\) is loss attributable to exact source/target identity;
- \(1-S_C\) is the fraction of within-endpoint relation information removed by
  the declared nuisance-control space.

The implementation exposes these quantities as
source_target_fe_retained_variance_fraction,
control_unique_variance_fraction_after_fe, and
total_unique_variance_fraction.

## Proposition 2 — nested-control survival

Let \(C_1\) be a baseline control set and let

\[
C_2 = [C_1, G]
\]

add one or more new controls \(G\) on the **same dyads**. Let the corresponding
residual focal relations be \(r_1\) and \(r_2\).

Define incremental survival as

\[
S_{2\mid 1}
=
\frac{\lVert r_2\rVert^2}{\lVert r_1\rVert^2}.
\]

Then

\[
S_{\mathrm{total},2}
=
S_{\mathrm{total},1} S_{2\mid 1}.
\]

For the climate case, \(G\) is continuous directed geographic
co-opportunity. Therefore

\[
1-S_{2\mid 1}
\]

is the exact fraction of previously unique C1 historical-relation information
removed by that additional geographic-control column space under the declared
linear projection.

This is a predictor-information statement. It does not imply that the removed
or retained component predicts a biological outcome.

## Proposition 3 — monotonicity under nested controls

For genuinely nested control spaces,

\[
\operatorname{col}(M_F C_1)
\subseteq
\operatorname{col}(M_F C_2).
\]

Orthogonal projection therefore gives

\[
0 \le \lVert r_2\rVert^2 \le \lVert r_1\rVert^2,
\]

and hence

\[
0 \le S_{2\mid 1} \le 1.
\]

Adding controls cannot increase residual focal-relation information except for
floating-point tolerance. The implementation checks this nested residual
sum-of-squares identity explicitly in the C1-to-C2 decomposition.

## Proposition 4 — scale invariance

For any nonzero scalar \(a\),

\[
S_F(ar)=S_F(r),
\qquad
S_C(ar)=S_C(r),
\qquad
S_{\mathrm{total}}(ar)=S_{\mathrm{total}}(r).
\]

Likewise, nonsingular rescaling of individual nuisance columns leaves their
column space unchanged and therefore leaves the residual focal relation
unchanged.

The information-survival fractions therefore describe geometry of the declared
predictor spaces rather than arbitrary measurement units.

## Signal concentration is a separate object

After the final declared control projection, define dyad-level signal mass

\[
w_{st}=r_{C,st}^{\,2}
\]

and normalized shares

\[
\pi_{st}
=
\frac{w_{st}}{\sum_{s,t} w_{st}}.
\]

TTF-Q reports the inverse-Herfindahl dyad concentration

\[
N_{\mathrm{dyad}}^{\mathrm{signal}}
=
\frac{1}{\sum_{s,t}\pi_{st}^{\,2}}.
\]

For source \(s\), aggregate

\[
\pi_s = \sum_t \pi_{st},
\]

giving

\[
N_{\mathrm{source}}^{\mathrm{signal}}
=
\frac{1}{\sum_s \pi_s^{\,2}},
\]

with an analogous definition for targets.

These values are bounded concentration summaries. They are **not effective
sample sizes**, do not determine standard errors, and do not substitute for
dependence-aware inference.

The known-truth benchmark deliberately holds total unique relation information
fixed while endpoint concentration changes. This demonstrates that information
amount and information concentration are distinct axes.

## Null qualification precedes detectability

Let \(A\) index a declared nuisance-heterogeneity regime and let \(p_0(A)\)
denote the synthetic rejection probability when the focal effect is
\(\beta=0\).

TTF-Q estimates that rejection probability by Monte Carlo simulation and
computes a Wilson interval

\[
[L_0(A), U_0(A)].
\]

Let \(\tau\) be the **predeclared Type-I qualification ceiling**. Define

\[
K(A)
=
\mathbf{1}\{U_0(A)\le\tau\}.
\]

For the current TTF-Q case study, the inferential test uses one-sided
\(\alpha_{\mathrm{test}}=0.025\), while the inherited design-qualification
ceiling is

\[
\tau=0.05.
\]

These are deliberately distinct quantities. Passing \(K(A)\) means that the
Monte Carlo Wilson upper bound satisfies the declared qualification ceiling;
it does **not** assert that the realized rejection probability equals exactly
0.025.

## Raw and evaluable minimum detectable effects

Let \(L_\beta(A)\) denote the Wilson lower bound on detection probability at
positive focal effect \(\beta\), and let \(\pi^\*\) be the target power.

The raw grid minimum detectable effect is

\[
D_{\mathrm{raw}}(A)
=
\min
\left\{
\beta>0:
L_\beta(A)\ge\pi^\*
\right\}.
\]

The **evaluable** grid MDE is then

\[
D_{\mathrm{eval}}(A)
=
\begin{cases}
D_{\mathrm{raw}}(A), & K(A)=1,\\
\text{undefined}, & K(A)=0.
\end{cases}
\]

This ordering is essential. A high rejection probability under a positive
effect is not interpretable as useful detectability when the same inferential
procedure violates its declared null-qualification ceiling in that nuisance
regime.

## Evaluable envelope

For a fixed dyadic predictor geometry, TTF-Q defines the response-blind
evaluable envelope as the joint object

\[
\mathcal{E}
=
\left(
S_F,\,
S_C,\,
S_{2\mid1},\ldots,\,
N_{\mathrm{dyad}}^{\mathrm{signal}},\,
N_{\mathrm{source}}^{\mathrm{signal}},\,
N_{\mathrm{target}}^{\mathrm{signal}},\,
\{K(A),D_{\mathrm{eval}}(A)\}_{A\in\mathcal{A}}
\right).
\]

Optional measurement-repeatability quantities may be appended when
response-blind resampling of the ecological relation is available.

The envelope intentionally does not collapse to one score. The known-truth
benchmarks show why:

- endpoint loss and nuisance redundancy can reduce identifiable information
  and worsen evaluable MDEs;
- concentrated signal can leave the total unique information fraction
  unchanged while violating the declared null-qualification ceiling;
- two different information-loss mechanisms can yield the same MDE vector;
- therefore information amount, information-loss location, signal
  concentration, null qualification, and conditional detectability answer
  different design questions.

## Known-truth validation

### Exact decomposition

The frozen decomposition benchmark prescribed \(S_F\), \(S_C\), and nested
geographic survival before results were generated. Across 243 scenarios, the
implementation recovered all prescribed fractions with maximum absolute error

\[
3.33\times 10^{-16}.
\]

Deliberately source- or target-concentrated signal was identified in every
predeclared comparison.

### Detectability and null qualification

The first frozen detectability benchmark revealed that raw MDEs can look
favorable even when the same estimator exceeds the declared Type-I ceiling.
That diagnostic was preserved rather than overwritten.

A second contract was then frozen before the calibration-repair results. It
retained the same scientific scenarios, increased null Monte Carlo precision,
and required \(K(A)=1\) before reporting an evaluable MDE.

Under that repaired definition:

- the broad reference design remained qualification-compatible, with evaluable
  MDE vector \(0.05/0.08/0.08/0.10\) across \(A=0/1/2/3\);
- endpoint loss and nuisance redundancy remained qualification-compatible but
  required larger effects;
- the deliberately source-concentrated design retained the same total unique
  relation fraction as the reference design but exceeded the declared
  \(\tau=0.05\) null ceiling at all tested nuisance amplitudes, so its
  evaluable MDE was undefined.

All four designs shared the same negative legacy-style
\((\beta=0.03,A=3)\) binary classification. The envelope separates the reason
for that common label.

## Interpretation boundary

None of these identities or synthetic-world results implies that the focal
ecological relation predicts a biological outcome.

TTF-Q characterizes the information geometry and the qualified inferential
domain available **before** that biological outcome is opened. The empirical
genetic-response firewall from relational v0.3 remains closed.
