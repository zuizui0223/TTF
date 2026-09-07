# TTF Method Specification v0.1

## 1. Scope

TTF (Transferable Turnover Fields) is a trait-agnostic framework for opportunistic georeferenced observations. Each record has a location \(x_{si}\) and trait value \(y_{si}\) for species \(s\). The trait only needs a defensible dissimilarity function \(d(y_i,y_j)\).

The inferential target is **cross-species transferability of within-species transition structure**, not concentration of all records in one pooled map.

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

To prevent high-amplitude or high-dimensional species from dominating solely through scale, standardize edge dissimilarities by rank **within species**:

\[
u_{se}=\frac{\operatorname{rank}_s(\delta_{se})-0.5}{|E_s|}.
\]

Ranks are the default measurement scale for the sharedness test. Raw trait amplitude remains a separate biological quantity and a calibration axis.

---

## 3. Layer 1 — opportunity-corrected boundary field

Let \(m_{se}\) be an edge midpoint and \(K_h\) a spatial kernel with bandwidth \(h\). Training species contribute equal total mass, so species with many observations do not dominate.

The continuous field is

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

The current implementation uses Euclidean coordinates. Geographic latitude/longitude should therefore be projected appropriately, or a future spherical graph/kernel implementation should be used.

---

## 4. Layer 2 — balanced repeats and descriptive recurrence

Repeated realizations are used to make the atlas stable to species and record inclusion.

The schedule must satisfy:

- no within-realization duplicate species;
- prospective species count per realization;
- long-run species inclusion counts differ by at most one;
- record-level balancing is required when repeated subsampling of records is used.

For fields \(B_r(x)\), a descriptive recurrence map is

\[
P_{\rm rep}(x)
=
P\{B_r(x)\text{ lies in the prespecified upper transition band}\}.
\]

**Important:** \(P_{\rm rep}\) is descriptive. It is not the cross-species sharedness estimand.

---

## 5. Layer 3 — structured permutation null

The primary null keeps fixed:

- species composition;
- observation coordinates;
- within-species graph nodes and edges;
- sampling schedule;
- environmental layers and predictor geometry;
- train/evaluation species split.

Only the trait-to-location assignment is broken, by shuffling trait rows **within species**.

For every permutation \(b\), the complete inferential path is rerun:

1. permute trait rows within species;
2. recompute edge dissimilarities;
3. rerank edges within species;
4. refit the training boundary field;
5. score held-out species.

A one-sided Monte Carlo p-value is

\[
p=\frac{1+\sum_b I(T^{(b)}\ge T_{\rm obs})}{B+1}.
\]

### Exchangeability caveat

Unrestricted within-species permutation is only valid when rows are exchangeable under the null. Repeated individuals, populations, years, or other dependence may require block-restricted permutation. TTF therefore accepts a per-record block label.

For pairwise genetic distances, do **not** shuffle distance-matrix cells. The appropriate future genetic interface should permute individual labels relative to coordinates while preserving the distance matrix's metric structure.

---

## 6. Layer 4 — primary sharedness estimand

Split species prospectively into disjoint training and evaluation sets.

Fit \(B^{(-S_{\rm eval})}\) using training species only.

For each evaluation edge, calculate mean boundary exposure along the edge:

\[
\hat b_{se}
=
\frac{1}{L_{se}}
\int_{e} B^{(-S_{\rm eval})}(x)\,d\ell .
\]

The implementation approximates this line integral with fixed midpoint quadrature points.

Within each held-out species,

\[
C_s=\operatorname{Spearman}(\hat b_{se},u_{se}).
\]

The primary statistic is the species-macro-average

\[
T=\frac{1}{|S_{\rm eval}|}\sum_{s\in S_{\rm eval}} C_s.
\]

This makes the scientific statement operational:

> boundary structure learned from other species predicts where this unseen species has its strongest within-species trait transitions.

No thresholded hotspot map is required for this test.

---

## 7. Layer 5 — predictor-space competition

Attribution is evaluated on the **same held-out species**, not by comparing in-sample fit.

A generic ridge-transfer implementation is provided for edge-level feature matrices. Intended prospective spaces include:

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

This supports statements such as “environmental information adds transfer skill beyond shared geography” more directly than winner-take-all model labels.

### Intrinsic IBD traits

For traits such as genetic distance, a within-species IBD expectation should be fit and cross-fit first, and sharedness should be tested on residual transition structure. That nuisance fit must be repeated under the null. This residualization layer is specified but not yet qualified in v0.1.

---

## 8. Layer 6 — mandatory sharedness × amplitude calibration

TTF is not considered interpretable until type-I error and power are separated across two axes:

- shared fraction \(\pi\);
- within-species signal amplitude \(A\).

The mandatory adversarial arm is

\[
\pi=0,\quad A>0.
\]

Every species may have a strong spatial boundary, but boundary locations are independent across species. Rejection must remain controlled as \(A\) increases.

The first implementation uses a circular geographic domain. Each species has a phase-shifted transition; private phases are uniform on the circle, while the prespecified fraction of shared species uses one common phase. This removes a privileged location from the \(\pi=0\) population while retaining strong spatial pattern.

Calibration must eventually be repeated on the actual empirical coordinate/graph geometry.

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

TTF explicitly separates three objects that can otherwise be confused:

1. **trait amplitude within species**;
2. **descriptive recurrence of a spatial field**;
3. **transferability across disjoint species**.

A method must not call (1) or (2) “sharedness” unless it also demonstrates (3).

The qualification suite therefore treats a strong-but-private transition world as a first-class null rather than as a nuisance discovered after the empirical analysis.
