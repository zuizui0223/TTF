# TTF Method Specification v0.2

## 1. Scope

TTF (Transferable Turnover Fields) is a trait-agnostic framework for opportunistic georeferenced observations. Each record has a location \(x_{si}\) and trait value \(y_{si}\) for species \(s\). The trait only needs a defensible dissimilarity function \(d(y_i,y_j)\).

The inferential target is **cross-species transferability of within-species transition structure**, not concentration of all records in one pooled map.

The core scientific statement is:

> transition structure learned from one set of species predicts where entirely unseen species exhibit strong within-species trait transitions.

---

## 2. Layer 0 — within-species transition geometry

For every species independently, construct a graph

\[
G_s=(V_s,E_s)
\]

using k-nearest neighbours, a distance threshold, or a prospectively fixed graph supplied by the analyst.

**Invariant:** an edge may never connect two species.

For each edge \(e=(i,j)\),

\[
\delta_{se}=d(y_{si},y_{sj}).
\]

Standardize edge dissimilarities by rank **within species**:

\[
u_{se}=\frac{\operatorname{rank}_s(\delta_{se})-0.5}{|E_s|}.
\]

Ranks are the default measurement scale for the sharedness test. Raw within-species trait amplitude is retained as a separate biological quantity and as an independent calibration axis.

---

## 3. Layer 1 — opportunity-corrected boundary field

Let \(m_{se}\) be an edge midpoint and \(K_h\) a spatial kernel with bandwidth \(h\). Training species contribute equal total mass, so species with many observations do not dominate.

\[
B(x)=
\frac{
\sum_{s\in S_{\rm train}} |E_s|^{-1}
\sum_{e\in E_s} K_h(x,m_{se})u_{se}+\lambda\mu_0
}{
\sum_{s\in S_{\rm train}} |E_s|^{-1}
\sum_{e\in E_s} K_h(x,m_{se})+\lambda
}.
\]

The denominator is geographic **edge opportunity** under the sampled graph. \(\lambda\) stabilizes low-opportunity locations and \(\mu_0=0.5\) is natural for rank scores.

For edge predictions, opportunity correction is applied at each quadrature point and only then averaged along the edge. The optimized fixed-geometry implementation is tested against the direct calculation so computational optimization cannot silently change the estimand.

The current implementation uses Euclidean coordinates. Latitude/longitude therefore needs an appropriate projection until a spherical graph/kernel interface is added.

---

## 4. Layer 2 — balanced repeats and descriptive recurrence

Repeated realizations can be used to make a descriptive atlas stable to species and record inclusion.

The schedule must satisfy:

- no within-realization duplicate species;
- prospective species count per realization;
- long-run species inclusion counts differ by at most one;
- record-level balancing when records are repeatedly subsampled.

For fields \(B_r(x)\), a descriptive recurrence map is

\[
P_{\rm rep}(x)=P\{B_r(x)\text{ lies in the prespecified upper transition band}\}.
\]

**Important:** \(P_{\rm rep}\) is descriptive. It is not the cross-species sharedness estimand and it is not a substitute for held-out-species inference.

---

## 5. Layer 3 — primary sharedness inference

### 5.1 Species-disjoint train/evaluation split

Split species prospectively into disjoint training and evaluation sets. Fit the boundary field using training species only.

For each evaluation edge, calculate mean boundary exposure along the edge:

\[
\hat b_{se}=\frac{1}{L_{se}}\int_e B^{(-S_{\rm eval})}(x)\,d\ell .
\]

Within each held-out species,

\[
C_s=\operatorname{Spearman}(\hat b_{se},u_{se}).
\]

The primary transfer statistic is

\[
T=\frac{1}{|S_{\rm eval}|}\sum_{s\in S_{\rm eval}}C_s.
\]

Thus sharedness is operationally defined as **out-of-species transferability**.

### 5.2 Conditional held-out-species null

The primary inferential null in v0.2 is

\[
H_0:E[C_s\mid B_{\rm train}]\le 0.
\]

The training field is frozen after learning from the training species. The evaluation species are then the inferential sampling units.

TTF uses a centered, studentized nonparametric bootstrap over the held-out species scores:

1. compute the observed held-out scores \(C_s\);
2. subtract their observed mean to impose a mean-zero null while preserving the empirical species-score distribution;
3. resample **species** with replacement;
4. compute the bootstrap mean and its studentized statistic;
5. compare the observed studentized mean to the centered bootstrap distribution.

No coordinates, traits, graphs, edge geometry, or within-species spatial pattern in evaluation species are altered.

The inferential claim is a **superpopulation claim about unseen species conditional on the trained field**. It is not a randomization test for the fixed finite set of evaluation species.

### 5.3 Why trait permutation is not the primary sharedness null

v0.1 used within-species trait-to-location permutation while holding coordinates and graph geometry fixed. That procedure remains useful for the distinct null of **trait-location exchangeability within species**.

However, prospective calibration showed that it is too narrow for the scientific null of no cross-species sharedness when every species may have a strong private spatial transition. Trait permutation destroys those private transitions. In the v0.1 circular-world pilot, the zero-shared rejection rate rose to 0.14 at the strongest private signal despite a prospective ceiling of 0.10, while positive-control power was 1.00.

Therefore v0.2 separates two questions:

- **primary sharedness test:** held-out-species conditional inference, preserving each species' spatial structure;
- **diagnostic trait-location test:** within-species trait permutation, subject to exchangeability assumptions.

This distinction is a methodological result, not an implementation detail.

### 5.4 Exchangeability caveat for diagnostic permutations

When diagnostic trait permutation is used, unrestricted within-species shuffling is valid only under row exchangeability. Repeated individuals, populations, years, or other dependence can require block-restricted permutation.

For pairwise genetic distances, do **not** shuffle distance-matrix cells. A genetic interface must preserve the distance matrix and permute individual labels relative to coordinates when that null is scientifically appropriate.

---

## 6. Layer 4 — qualification by sharedness × amplitude

TTF separates two axes that pooled concentration statistics can confound:

- shared fraction \(\pi\);
- within-species spatial signal amplitude \(A\).

The mandatory adversarial null arm is

\[
\pi=0,\quad A>0.
\]

Every species may have a strong spatial transition, but transition locations are private. Type-I error must remain controlled as \(A\) increases.

The first synthetic implementation uses a circular geographic domain. Private species have uniformly random boundary phase. A prespecified fraction of species shares a common phase. This orthogonalizes transition strength from cross-species sharing.

### v0.1 trait-permutation pilot

50 worlds per cell, 199 permutations, 40 species, 60 records/species:

- maximum zero-shared positive-amplitude rejection: **0.14** — fail;
- full-sharing moderate-amplitude power: **1.00** — pass.

This isolated a null-construction failure rather than a failure of the edge-level transfer estimand.

### v0.2 held-out-species pilot

100 worlds per cell, 1,999 species-bootstrap resamples, same 40-species/60-record geometry:

- zero-shared rejection at \(A=0.5,1,2,3\): **0.01, 0.07, 0.05, 0.07**;
- maximum zero-shared positive-amplitude rejection: **0.07** — passes provisional ceiling 0.10;
- full-sharing, moderate \(A=2\) power: **0.96** — passes floor 0.80;
- full-sharing \(A=3\) power: **1.00**.

This is a **provisional idealized-world qualification**, not publication-grade calibration. With 100 worlds, Wilson confidence bounds on a 0.07 rejection rate remain wider than the desired final precision.

The next gates are higher-replicate type-I calibration and semi-synthetic calibration on actual opportunistic sampling geometry.

---

## 7. Layer 5 — predictor-space competition

Attribution is evaluated on the **same held-out species**, not by comparing in-sample fit.

Intended prospective spaces include:

- \(D\): smooth geographic distance / intrinsic IBD expectation;
- \(G\): geographic barriers or spatial structure;
- \(E\): environmental contrast / threshold structure;
- \(GE\): geography + environment;
- \(GER\): geography + environment + resistance.

Primary attribution should use paired held-out increments such as

\[
\Delta_{E|G}=T_{GE}-T_G,\quad
\Delta_{G|E}=T_{GE}-T_E,\quad
\Delta_{R|GE}=T_{GER}-T_{GE}.
\]

The same species split must be used for every competing predictor space.

### Intrinsic IBD traits

For traits such as genetic distance, first estimate a within-species IBD expectation and test transferable structure in residuals. The nuisance model must itself be cross-fit so an evaluation species is not used to construct the predictor against which it is scored.

This residualization layer is not yet qualified.

---

## 8. Layer 6 — actual-geometry calibration

Idealized circular-world qualification is necessary but not sufficient. Before empirical deployment, repeat semi-synthetic calibration using the intended empirical:

- species counts;
- record counts;
- coordinates;
- within-species graph geometry;
- opportunity structure;
- train/evaluation constraints.

Trait values remain synthetic. This establishes the dataset-specific detection floor without letting the observed trait outcome determine the calibration world.

Report a detection surface or minimum detectable amplitude as a function of shared fraction rather than reducing identifiability to one power number.

---

## 9. Layer 7 — observation-bias controls

Sampling availability should be frozen before outcome analysis. Candidate controls include:

- observation density / platform activity;
- licensing availability;
- annotation or measurement availability;
- observer caps;
- deletion of highest-effort cells;
- location-blind trait measurement;
- negative-control availability surfaces.

Even after these controls, the claim ceiling is **reproducible transition structure within the measurable sampling frame**, not proof that opportunistic sampling is unbiased.

---

## 10. Design decisions inherited from the RGFCA lessons

TTF explicitly separates four objects that can otherwise be confused:

1. **trait amplitude within species**;
2. **descriptive recurrence of a spatial field**;
3. **trait-location exchangeability within species**;
4. **transferability across disjoint species**.

Only (4) is the primary sharedness target.

The v0.1 failure is therefore retained as part of the method's audit trail: preserving sampling geometry during trait permutation was not enough, because the scientific no-sharedness null also needs to allow strong private within-species spatial structure.

A method is not qualified merely because it controls type I in a no-structure world. It must also control type I when every species has strong structure in different places and must recover a prospectively defined shared positive control.
