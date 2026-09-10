# TTF-M / TTF-C — Transferable Mismatch and Coupling Fields

Status: **method development with idealized direction-free qualification passed; Q4 directionality v0.1 failed and is under a new prospective strict-private qualification design**. This extension begins after the frozen TTF v0.1–v0.12 history and does not reinterpret any earlier failure or PASS. It opens no empirical flower-colour outcome.

## 1. Purpose

Core TTF asks whether a field of strong within-system state turnover learned from training systems predicts turnover in entirely unseen systems.

TTF-M generalizes the state being transferred from one field `Y_s(x)` to a declared relationship between two co-located fields `A_s(x)` and `B_s(x)`.

The extension separates two questions that must not be conflated:

1. **Mismatch-magnitude transition** — where does the amount of mismatch between A and B change?
2. **Coupling-relation transition** — where does the A–B relationship itself change, even when mismatch magnitude is unchanged?

The primary scientific principle is:

> A transition structure is system-general only to the extent that a turnover-intensity field learned without the evaluation systems predicts within-system turnover in those unseen systems.

This is deliberately broader than “the same exact boundary occurs in every system.” TTF can detect a non-uniform **distribution of transition locations** when that distribution induces a transferable turnover-intensity field. Exact-front sharing is a special case.

A useful semantic hierarchy is therefore

```text
exact deterministic front sharing
    subset of
shared front-density / turnover-intensity structure
    subset of
transferable transition structure detected by TTF
```

Independent exact front locations are not automatically a TTF null. If those locations are drawn from a common non-uniform spatial density, their expected turnover field can remain predictable in held-out systems.

## 2. Paired observation support

For system `s`, observations must be paired on the same predeclared support:

\[
\{(x_{si}, A_{si}, B_{si})\}_{i=1}^{n_s}.
\]

Rows may not be independently re-matched after inspecting mismatch or transfer outcomes. Missingness, interpolation, nearest-neighbour matching, temporal windows, and aggregation must be declared before qualification because they change the paired estimand.

`PairedSpeciesSample` is the current implementation carrier. The word `species` is inherited from core TTF; a future unit-general API may rename it to `system` without changing the estimand.

## 3. TTF-M: mismatch-magnitude field

Declare a pointwise mismatch function

\[
M_s(x)=m\{A_s(x),B_s(x)\}\ge 0.
\]

Examples include absolute phenological lag, phenotype–environment distance, consumer–resource trait distance, or a predeclared incompatibility score.

For a fixed within-system graph edge `e=(i,j)`, define raw mismatch turnover

\[
q^M_{se}=d_M\{M_{si},M_{sj}\},
\]

then apply the same within-system rank standardization as core TTF:

\[
u^M_{se}=\operatorname{rank}_s(q^M_{se}).
\]

The ordinary TTF field learner and held-out transfer statistic are then applied to `u^M` without modification.

### Interpretation

A positive TTF-M transfer result means that the mismatch-turnover intensity learned in training systems predicts mismatch turnover in unseen systems.

Because core TTF uses undirected edges, TTF-M detects a **mismatch-transition field**, not by itself the direction `low -> high mismatch`. Calling a transition a mismatch *rise* or *breakdown* requires a separately frozen orientation rule.

## 4. TTF-C: coupling-relation field

Mismatch magnitude can be too coarse. Let a declared relation representation be

\[
R_s(x)=r\{A_s(x),B_s(x)\}.
\]

When A and B live in a common vector space, the default is

\[
R_s(x)=A_s(x)-B_s(x).
\]

Coupling turnover on edge `e=(i,j)` is

\[
q^C_{se}=d_R\{R_{si},R_{sj}\},
\]

with within-system rank score

\[
u^C_{se}=\operatorname{rank}_s(q^C_{se}).
\]

Again, the existing held-out TTF field learner is used unchanged.

### Why TTF-C is stricter than mismatch magnitude

If A and B undergo the same common shift,

\[
\Delta A = \Delta B,
\]

then `R=A-B` is unchanged and TTF-C has no coupling transition even though both component fields may show strong turnover.

Conversely, equal mismatch magnitudes can conceal a changed relation. In two dimensions,

\[
R_1=(1,0), \qquad R_2=(0,1)
\]

have equal norms but different coupling states. TTF-M based only on `||R||` is constant, whereas TTF-C detects the relational transition.

Thus TTF-C is not a residualized version of TTF-M. It changes the target representation from mismatch magnitude to the declared coupling relation itself.

## 5. What counts as "breakdown"

TTF-C measures **coupling change**, not automatically deterioration.

A coupling-change field may represent decoupling, recoupling, a rotation or qualitative relation change, asymmetric movement of A and B, or opposite movement of A and B.

To claim specifically a *breakdown*, an application must additionally predeclare how the transition is oriented by coupling quality, mismatch level, compatibility, fitness, or another biologically justified criterion.

Q4 formalizes one such design with a signed orientation coordinate `g_s(x)` fixed independently of observed mismatch outcomes. Its v0.1 qualification failed and remains failed; see `Q4_DIRECTIONALITY_FAILURE_DIAGNOSIS.md`.

## 6. Difference from SDM and ordinary mismatch regression

### SDM / state model

Typical target:

\[
Y_s(x) \sim X(x).
\]

Question: where is a state, occurrence, abundance, trait value, or suitability high?

### Ordinary mismatch model

Typical target:

\[
M_s(x) \sim X(x).
\]

Question: which predictors explain mismatch magnitude?

### TTF-M

Target:

\[
\text{held-out transfer of } \Delta M_s(e).
\]

Question: does a mismatch-turnover field learned in some systems predict mismatch turnover in unseen systems?

### TTF-C

Target:

\[
\text{held-out transfer of } \Delta R_s(e).
\]

Question: does an A–B relation-turnover field learned in some systems predict relation turnover in unseen systems?

TTF-M/TTF-C are therefore not replacements for SDMs. They address cross-system reproducibility of transition geometry/intensity rather than the response surface itself.

## 7. Qualification families and current state

### Q0 — coupled shared component transition, no mismatch/coupling transition

A and B share the same spatial transition in every system. Both component fields have a highly transferable boundary, but `A-B` and mismatch magnitude remain constant.

Required behavior: TTF-M and TTF-C type-I remain controlled. **Passed** in the frozen idealized qualification.

### Q1 — strict-private mismatch/relation transitions

A valid TTF private null must contain strong system-specific transitions while making the population expected turnover-intensity field spatially constant on the transfer domain. The idealized circular full-phase construction has this property by rotational symmetry.

Required behavior: held-out TTF-M and TTF-C type-I remain controlled despite strong within-system structure. **Passed** in the frozen idealized qualification.

A finite-line construction with front centers restricted to a common central band is **not** equivalent to this strict-private null; it creates a shared front-density zone.

### Q2 — shared mismatch-magnitude transition field

A and B decouple at a shared spatial transition so that mismatch magnitude changes across systems.

Required behavior: TTF-M has predeclared power while preserving Q0/Q1 validity. **Passed** in the frozen idealized qualification.

### Q3 — shared relational rotation at constant mismatch magnitude

The relative state `R=A-B` changes orientation across a common boundary while `||R||` remains constant.

Required behavior: TTF-M stays null while TTF-C detects relation turnover. **Passed**: the frozen gate produced TTF-M 0/500 versus TTF-C 500/500.

### Q4 — directionality: breakdown versus recoupling

Q4 requires both a transferable TTF-C relation field and a held-out oriented mismatch sign. Q4 v0.1 **failed formally** because its declared PRIVATE cell exceeded the type-I ceiling.

Post-outcome diagnosis showed that this PRIVATE cell drew independent exact front locations only from a central subdomain. It therefore contained a non-uniform transferable front-density field and was not a clean null for TTF's actual estimand. The failure remains immutable; a new strict-private protocol is required rather than threshold retuning.

### Q5 — density and geometry stress

Repeat the required estimand/direction controls on heterogeneous opportunistic geometries and density-scaled graph rules, preserving the lessons of TTF v0.8–v0.11. Q5 remains blocked until Q4 directionality semantics are resolved and prospectively requalified.

## 8. Primary endpoint and inference status

For each direction-free estimand, retain the core TTF primary endpoint:

\[
T=|S_{eval}|^{-1}\sum_s
\operatorname{Spearman}(\hat b_{se},u_{se}).
\]

Use completely system-disjoint training/evaluation sets. Null construction must preserve strong system-private mismatch/coupling structure without creating a non-constant population turnover-intensity field unless that field is intentionally being tested as an alternative.

The frozen idealized qualification used 40 systems × 60 records, 20/20 train/evaluation, density-scaled `k=9`, 500 worlds per cell, 1,999 held-out-system bootstrap resamples, alpha 0.05, Wilson 95% type-I upper ceiling 0.10, and power lower floor 0.80. TTF-M and TTF-C both passed that declared family.

Q4 directionality has a separate failed formal result and is not licensed by the idealized direction-free PASS.

## 9. Application classes

Potential applications include plant flowering vs pollinator activity phenological mismatch; floral phenotype vs local pollinator functional composition; phenotype vs environmental optimum; consumer vs resource trait matching; realized occurrence vs independently estimated suitability; host–symbiont or mutualist compatibility; genotype/phenotype vs local selective environment; and spatial or spatiotemporal coupling fronts under climate change.

These are examples, not validated empirical use cases.

## 10. Claim ceiling

The current branch supports the following bounded statements only:

- TTF-M and TTF-C have distinct implemented estimands;
- their frozen balanced idealized qualification passed its declared type-I/power criteria;
- TTF-C detects relation change that can be invisible to mismatch magnitude;
- Q4 v0.1 directionality qualification failed and does not license breakdown/recoupling use;
- exact-front independence is not sufficient to define a TTF null when front-location density itself creates transferable turnover intensity.

It does **not** establish empirical validity in any ecological system, causal mismatch mechanisms, directional breakdown qualification, heterogeneous-geometry robustness, fitness consequences, or any flower-colour/pollination conclusion. It also does not imply that `A-B` is appropriate for heterogeneous state spaces without a justified adapter, that a transferable transition field has one shared environmental driver, or that TTF-M/TTF-C supersede SDMs, mismatch regressions, joint species models, or causal interaction models.

Frozen TTF v0.1–v0.12 results remain unchanged and are not evidence for this extension.