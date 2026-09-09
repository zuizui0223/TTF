# TTF-M / TTF-C — Transferable Mismatch and Coupling Fields

Status: **prospective method-development specification**. This extension begins after the frozen TTF v0.1–v0.12 history and does not reinterpret any earlier failure or PASS. It opens no empirical flower-colour outcome.

## 1. Purpose

Core TTF asks whether locations of strong within-system state turnover learned from training systems predict turnover in entirely unseen systems.

TTF-M generalizes the state being transferred from one field `Y_s(x)` to a declared relationship between two co-located fields `A_s(x)` and `B_s(x)`.

The extension separates two questions that must not be conflated:

1. **Mismatch-magnitude transition** — where does the amount of mismatch between A and B change?
2. **Coupling-relation transition** — where does the A–B relationship itself change, even when mismatch magnitude is unchanged?

The primary scientific principle remains unchanged:

> A transition is system-general only to the extent that a field learned without the evaluation systems predicts within-system transition in those unseen systems.

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

The ordinary TTF boundary learner and held-out transfer statistic are then applied to `u^M` without modification.

### Interpretation

A positive TTF-M transfer result means that locations where mismatch magnitude changes in training systems predict mismatch-transition locations in unseen systems.

Because core TTF uses undirected edges, TTF-M detects a **mismatch-transition front**, not by itself the direction `low -> high mismatch`. Calling a front a mismatch *rise* or *breakdown* requires a separately frozen side-label or orientation rule.

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

A coupling-change boundary may represent:

- decoupling;
- recoupling;
- a rotation or qualitative change in relation;
- asymmetric movement of A and B;
- opposite movement of A and B.

To claim specifically a *breakdown front*, an application must additionally predeclare how the sides of the front are ordered by coupling quality, mismatch level, compatibility, fitness, or another biologically justified criterion.

This prevents direction-free turnover from being overinterpreted as deterioration.

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

Question: do mismatch-transition locations learned in some systems predict where mismatch changes in unseen systems?

### TTF-C

Target:

\[
\text{held-out transfer of } \Delta R_s(e).
\]

Question: do locations where the A–B relation changes transfer to unseen systems?

TTF-M/TTF-C are therefore not replacements for SDMs. They address cross-system reproducibility of transition geometry rather than the response surface itself.

## 7. Prospective synthetic qualification families

No empirical application should be opened before the following families are implemented and frozen.

### Q0 — coupled shared component transition, no mismatch/coupling transition

A and B share the same spatial transition in every system. Both component fields have a highly transferable boundary, but `A-B` and mismatch magnitude remain constant.

Required behavior:

- core component TTF may be positive;
- TTF-M type-I must remain controlled;
- TTF-C type-I must remain controlled.

This is the critical protection against calling ordinary shared turnover "decoupling".

### Q1 — private mismatch transitions, no shared front

Each system has a strong mismatch/coupling transition at an independently sampled location.

Required behavior: held-out TTF-M and TTF-C type-I remain controlled despite strong within-system structure.

### Q2 — shared mismatch-magnitude front

A and B decouple at a common boundary so that mismatch magnitude changes across systems.

Required behavior: TTF-M has predeclared power while preserving Q0/Q1 validity.

### Q3 — shared relational rotation at constant mismatch magnitude

The relative state `R=A-B` changes orientation across a common boundary while `||R||` remains constant.

Required behavior:

- TTF-M remains null or low-power by construction;
- TTF-C detects the shared coupling transition.

This is the discriminating qualification for the added TTF-C estimand.

### Q4 — shared recoupling / reverse-side control

The same coupling-change geometry occurs but the high-quality/low-mismatch side is reversed relative to a breakdown scenario.

Required behavior: direction-free TTF-C may detect the boundary, but no generic "breakdown" claim is licensed without the separate side-label rule.

### Q5 — density and geometry stress

Repeat Q0–Q4 on heterogeneous opportunistic geometries and density-scaled graph rules, preserving the lessons of TTF v0.8–v0.11. A PASS on idealized regular coordinates cannot substitute for this stress layer.

## 8. Planned primary endpoints

For each estimand, retain the core TTF primary endpoint:

\[
T=|S_{eval}|^{-1}\sum_s
\operatorname{Spearman}(\hat b_{se},u_{se}).
\]

Use completely system-disjoint training/evaluation sets. Null construction must preserve strong system-private mismatch/coupling structure rather than erase it by naive pointwise permutation.

Numerical type-I/power thresholds are **not set by this specification**. They must be frozen in a later qualification protocol before synthetic outcomes are opened.

## 9. Application classes

Potential applications include:

- plant flowering vs pollinator activity phenological mismatch;
- floral phenotype vs local pollinator functional composition;
- phenotype vs environmental optimum;
- consumer vs resource trait matching;
- realized occurrence vs independently estimated suitability;
- host–symbiont or mutualist compatibility;
- genotype/phenotype vs local selective environment;
- spatial or spatiotemporal coupling fronts under climate change.

These are examples, not validated use cases.

## 10. Claim ceiling

TTF-M/TTF-C currently establish only software-level estimand definitions and prospective qualification targets.

They do not yet establish:

- calibrated type-I or power for mismatch/coupling transfer;
- empirical validity in any ecological system;
- that a coupling transition is causal;
- that a detected coupling-change front is deterioration rather than recoupling;
- that `A-B` is appropriate for heterogeneous state spaces without a justified adapter;
- that a shared transition implies one shared environmental driver;
- that TTF-M/TTF-C supersede SDMs, mismatch regressions, joint species models, or causal interaction models.

Frozen TTF v0.1–v0.12 results remain unchanged and are not evidence for this extension.