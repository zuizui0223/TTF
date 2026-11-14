# Unknown-variance historical-signal design benchmark

> Implementation workflow: superpowers:executing-plans; deterministic test-first mathematics. User has delegated continuation; no empirical response opening.

**Goal:** Determine whether the existing source geometry retains historical information after absorbing arbitrary species and study intercepts, without needing their variance components or a known residual variance.
**Architecture:** A source-only projection eliminates [21 controls, intercept, species indicators, study indicators]. A hypothetical Gaussian homoscedastic independent residual has unknown scale and yields an exact Student t pivot. No original rates are accepted by the implementation. Reuse byte-locked input reconstruction; never edit frozen source modules.
**Stack:** NumPy, SciPy; existing pinned Ubuntu24/Python3.12.15 environment on GitHub Actions. Local deterministic tests may use installed numerical versions; primary original-source run uses pinned versions.
**Spec:** This file.

## Approaches and decision

- Estimate mixed-model variance components: realistic future candidate, but estimation uncertainty requires a new empirical or simulation programme; not chosen here.
- Known covariance-shape with unknown overall scale: removes only one unknown scale and retains unverified species/study variance ratios; not sufficient for this task.
- **Selected:** absorb species and study intercepts, allow unknown iid residual scale. Gives a finite-sample conditional t reference under Gaussian independent homoscedastic residuals, at the cost of within-study/species information and degrees of freedom. This is a design benchmark, not a replacement PASS for Haldane v1.

## Fixed scientific choices

- Source X/H: the exact original 165x21 and 165 arrays; original 91/74 row split, 29/27 species, 48/45 studies.
- Test-side geometry only: 74 rows, 27 species, 45 studies. Training species-equal means/SDs scale h and controls, unchanged from the prior oracle benchmark.
- W = [intercept, standardized controls, full species indicators, full study indicators]. Column-normalized SVD, relative rank cutoff 1e-10, fixed before source evaluation.
- r=(I-P_W)h. df=74-rank(W)-1. Stop diagnostic as nonidentifiable if r is numerically zero; no t comparison if df<=0.
- Hypothetical y=Wb+gamma*h+epsilon, epsilon~N(0,sigma^2 I), sigma unknown. Group effects are arbitrary unknown additive intercepts, not assumed to have known variances or distributions.
- u=r/||r||, M=I-P_W-uu'. Under gamma=0, (u'y)/sqrt(y'My/df) ~ t_df. Group effects cancel; residual independence, homoscedasticity, Gaussianity and correct linear mean remain assumptions.
- Nominal two-sided alpha=1/60. Power is noncentral-t with ncp=gamma*||r||/sigma. Report gamma=0,.25,.5,1 and 80% detectable gamma. For comparison only, set hypothetical sigma=.5 from the old iid noise component; it is NOT estimated biological noise.
- Report information retention relative to controlling the original 21 variables alone; degrees of freedom; projection leakage; residual-signal concentration by species/study; no per-site weights/arrays.
- No actual y interface, no random number draws, no new null worlds, no outcome opening, no old thresholds/hash/source edits.

## Tasks and verification

- [x] Core projection and noncentral-t power: deterministic red test, minimal implementation, green tests; invariance under arbitrary group offsets, covariance of group effects, global scale, duplicated nuisance columns; fail finite/shape/rank guards.
- [x] Source runner: exact original code and raw climate/X/H identity; unchanged terminal receipt; additive new files only. Test source hash guards and no random/empirical-opening calls.
- [ ] Pinned workflow: run new deterministic tests plus existing oracle/analytic/numeric tests; one original-source reconstruction. Export aggregate metrics only.
- [ ] Read artifact, verify digest, record interpretation and actual CI status. Do not label design validity as real biological validation.

Local deterministic tests: 17 passed after two separate failing implementation-existence test runs. Local NumPy 2.3.5/SciPy 1.17.0; primary source run retains original pinned versions. Direct local git clone was unavailable due to network DNS; connected GitHub API used for source reads and atomic branch update, with no main changes.

## Review focus

1. Absorption can eliminate the history signal or residual degrees of freedom: report nonidentifiable, do not retune rank or drop controls.
2. Unknown group variance cancellation is NOT robustness to arbitrary residual covariance, random slopes, spatial/phylogenetic dependence or heteroscedasticity.
3. Same pair rows may not be independent residuals in actual experiments; this diagnostic does not establish that assumption.
4. A common coefficient and conditional association are not unknown-species predictive gain or causal adaptation/plasticity/speciation.
5. Entire repository CI is separate from dedicated deterministic tests; report pending/failures rather than claiming whole-repository success.

References: R stats summary.lm official manual (residual sigma and coefficient t); SciPy stats.t and stats.nct official documentation (finite-sample and noncentral-t probabilities). The group projection proof follows directly from orthogonal projection and independent Gaussian linear/quadratic forms.
