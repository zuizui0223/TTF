# Prospective qualification and abstention for cross-species transfer claims

**Methods manuscript v0.1 — methods-only refactor**

Target journal: *Methods in Ecology and Evolution*

## One-sentence contribution

Cross-species transfer should be treated as a claim that must first be shown to be evaluable on the exact deployment geometry; TTF provides a prospective workflow that can authorize one empirical opening or terminate as `NOT_EVALUABLE` without outcome-dependent rescue.

## Abstract

Ecologists increasingly use models trained in one set of species, populations, sites or systems to make claims about structure shared beyond the training data. Standard held-out prediction is necessary for such claims but is not sufficient when the deployment geometry itself can make a statistic anti-conservative, weakly identifiable or observationally unsupported. We introduce TTF, a prospective qualification-and-abstention framework for cross-unit transfer. TTF separates method validity, deployment-geometry detectability, observational admissibility and empirical effect estimation into sequential gates. Candidate estimands, nuisance controls, synthetic worlds, calibration thresholds and stopping rules are frozen before held-out outcomes are opened. The exact empirical geometry must then pass prospective Type-I and power qualification; designs that fail terminate as `NOT_EVALUABLE` rather than being retuned against the observed response. We demonstrate the framework on cross-species spatial genetic differentiation. A development geometry failed the predeclared cross-species error-control margin despite passing a separate within-species detectability test, correctly preventing an empirical opening. A fresh phylogatR panel then passed response-blind intake and character-support filters, leaving 211 species with an inherited 103/108 training/evaluation split. Its exact geometry passed the frozen cross-species synthetic gate (maximum private-null Wilson 95% upper bound 0.0681; shared moderate-signal Wilson 95% lower bound 0.9714), authorizing one empirical analysis. Post-IBD transfer to unseen species was not detected (T = 0.0326, p = 0.7393), while a separately qualified within-species self-detectability statistic exceeded its frozen null reference (p = 0.0090). The substantive genetic result is secondary to the methodological demonstration: a transfer workflow can distinguish a valid negative empirical result from a design that was never capable of supporting the claim. TTF turns abstention from an informal caveat into an executable endpoint.

## 1. The problem

A held-out test answers whether a fitted model predicts unseen observations under a chosen split. It does not by itself establish that the intended scientific transfer claim was identifiable on the realized geometry.

Three failure modes are especially important in ecological transfer problems:

1. **Private structure can mimic shared structure.** Strong within-unit autocorrelation may generate apparent transfer under a null that destroys the dependence structure that should have been preserved.
2. **Deployment geometry can change calibration and power.** A procedure that works under one record density, graph scale or spatial support can fail under another even when the estimator is unchanged.
3. **Observational support can fail after method qualification.** The inferential machinery may be valid in principle while the realized empirical data do not support the predeclared estimand.

A valid framework therefore needs a state between `SIGNIFICANT` and `NON_SIGNIFICANT`: the analysis must be able to say that the requested claim was not evaluable under the frozen design.

## 2. TTF design

TTF enforces the sequence

```text
estimand freeze
  -> response-blind deployment geometry
  -> observational/admissibility gate
  -> exact-geometry synthetic qualification
  -> one authorized empirical opening
  -> terminal interpretation
```

Each transition has an explicit machine-readable receipt. Failure at a gate prevents later response access.

### 2.1 Species-disjoint transfer

The unit of transfer is the biological unit, not the individual edge, record or pairwise comparison. Training and evaluation species are disjoint. For evaluation species (s),

```text
C_s = association(prediction learned without s, observed structure in s)
T   = mean_s C_s.
```

The particular response can change across applications; the held-out-unit logic does not.

### 2.2 Nuisance-safe response construction

For pairwise genetic data, target edges share endpoint localities with other edges. TTF therefore fits the within-species isolation-by-distance nuisance relation for target edge `(i,j)` using only edges that share neither `i` nor `j`. The primary response is the residual spatial structure after this endpoint-safe nuisance fit.

### 2.3 Exact-geometry qualification

Synthetic qualification is run on the exact empirical graph geometry that would be used after opening the response. The private-null family preserves within-species spatial structure but contains no transferable shared field. Shared-signal worlds add a known transferable component.

A deployment is authorized only when both error control and minimum power pass their prospectively frozen confidence-bound criteria.

### 2.4 Abstention

If calibration, power, character support or observational support fails, TTF returns a named `NOT_EVALUABLE` state and blocks the empirical statistic. A failed frozen family is not repaired by trying a new bandwidth, graph, threshold, null or response transformation after seeing the result.

### 2.5 Finite candidate families

When an estimand requires a candidate mechanism, transport or relational predictor, the candidate family is declared finite before the confirmatory response is inspected. Exhausting that family closes the program. A successor must be a newly named hypothesis with a new development family rather than an additional member appended after failure.

## 3. Demonstrations

### 3.1 Development geometry: a useful failure

The development cross-species Gate-D failed its frozen private-structure error-control bound: the Wilson upper bound was 0.1003347533 against a required maximum of 0.10. The shared-signal power requirement passed.

A separate within-species self-detectability qualification passed on the declared synthetic family.

The correct endpoint was therefore not a biological null. The development dataset was **not evaluable for the cross-species empirical claim under that design**, and its genetic outcomes remained unopened.

This near-threshold failure is important because it shows that TTF does not convert a nearly acceptable calibration into a usable empirical test.

### 3.2 Observational admissibility can fail after synthetic qualification

In an earlier flower-colour deployment, the method-level synthetic qualification passed, but the independently frozen measurement pipeline left only 189 of 250 species above the required support threshold and no authorized empirical geometry was retained. Trait values, pairwise colour distances and the empirical TTF statistic remained unopened.

This is a different failure mode from poor Type-I calibration. It demonstrates why observational admissibility must be separated from method validity.

### 3.3 Fresh genetic deployment

A fresh phylogatR COI/COX1 intake froze 250 species using response-blind geometry. Character-support filtering left 211 survivors and inherited the original deterministic split, yielding 103 training and 108 evaluation species.

The exact survivor geometry passed the prospective synthetic qualification:

- maximum private-null Wilson 95% upper rejection bound: **0.0681**;
- shared moderate-signal Wilson 95% lower power bound: **0.9714**.

A separate within-species self-detectability gate also passed, authorizing one empirical response opening.

The primary post-IBD cross-species transfer statistic was **T = 0.0326, p = 0.7393**. The separately calibrated within-species self statistic had **p = 0.0090** relative to its frozen structural null.

The methodological interpretation is not that TTF “found no pattern.” The workflow established that this was an **authorized negative empirical result**, rather than a negative result produced by an unqualified geometry.

## 4. What the framework changes

TTF changes the interpretation of negative and failed analyses.

| State | Meaning |
| --- | --- |
| Qualification FAIL | The requested empirical transfer claim is not licensed under the design. |
| Observational gate FAIL | The method may be valid, but the realized data cannot support the predeclared estimand. |
| Qualification PASS + empirical non-detection | A valid negative result within the qualified domain. |
| Qualification PASS + empirical detection | Evidence for transfer under the frozen estimand and null family. |

This distinction is the primary methodological output.

## 5. Scope and limitations

Qualification is conditional on the declared synthetic world family. Passing does not prove that every biologically plausible source of private structure has been represented. TTF therefore requires the simulation family and its omissions to be explicit.

Likewise, abstention protects the confirmatory interpretation but does not identify the biological reason a geometry or observational pipeline failed. Mechanistic follow-up belongs to a newly declared development program.

The framework is most useful when the scientific claim is explicitly about reuse across biological units and when dependence, geometry or measurement support can invalidate a naive held-out test.

## 6. Reproducibility contract

Every confirmatory claim should be traceable through

```text
frozen rule
  -> source identity
  -> response-blind intake
  -> qualification implementation
  -> immutable gate receipt
  -> authorized empirical result
  -> manuscript claim
```

Failed gates and closure receipts are retained as first-class outputs because they establish that the workflow can actually stop.

## 7. Planned repository refactor

The butterfly specialization ecology paper is being migrated to `zuizui0223/chocho`. Once that repository independently passes its paper-specific tests, butterfly/Lepidoptera ecology scripts, exploratory receipts and submission materials will be removed from the TTF methods branch. Frozen methodological history needed to document qualification or abstention will remain only when it contributes directly to the methods argument.
