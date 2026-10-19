# LGM resource-history specificity: exploratory v0.5 ecological verdict (19 Oct 2026)

**Status: climate-analogue predictor screening only. Negative for known-host identity advantage under a relevant botanically sourced control. No genetic or historical occupancy confirmation.**

## Core evolutionary ecology question

Do insect species living in the same glacial landscape develop different genetic breaks because each depended on its own larval host plants, thereby inheriting different historical resource connectivity?

The current proxy is **not yet actual host connectivity**. Each insect's source-hashed modern sampling geometry is evaluated at five interior positions along paired locations, using CHELSA-TraCE21k 0 BP vs LGM (21 ka; 4 climate variables) distances to modern observational host or insect climate reference sets.

## Fixed input layers and successive results

| Screening step | Source | Result | What it does/does not show |
|---|---|---:|---|
| 397 climate-finite species, modern host references independently separated from focal sample locations and matched to insect reference count | v0.3, draft PR #147 | Median host-history unique information fraction **3.714%**, conditional on own-history, both present climate mismatches, and geographic distance | Some historical climatic predictor variation remains. **Not genetic R² or glacial host occupancy** |
| Independent same-insect-niche reference split, 347 climate-finite species | v0.4, draft PR #148 | Actual known plant hosts > pseudo-host built from another insect sample in **281/347** species | Actual plant climate differs from insect climate; this alone does **not** demonstrate larval-host specificity |
| Modern-climate-matched other real plants, 335 source-qualified species / 333 climate-finite | **v0.5, draft PR #149** | Actual host > alternative plant in **147/333**; actual raw residual variance > alternative in **142/333** | **No overall host-identity advantage** under this alternative botanical background |

The negative v0.5 result was computed from **eight source-frozen alternative plant taxon groups/species**, all selected using **0 BP climate signatures before examining LGM scores**, excluding the recorded focal host IDs and genera. An unlisted plant is **not certified to be an unused larval host**.

## Correct paired interpretation of the v0.5 numbers

- Separate medians: actual-host information fraction **3.9391365%**, alternative-plant **3.8016618%**.
- **The actual paired species-level contrast is not the subtraction of those marginal medians.** Median `actual - alternative`: **−0.1685039 percentage points**. Thus 186/333 favor the alternative rather than the known host.
- Median paired raw unexplained climate-history variance is also lower for actual hosts (**−0.000351114**), preventing the interpretation that a positive host effect was hidden solely by distinct fraction denominators.
- Same 333 insects: v0.4 insect-sample pseudo-host comparison gave **269/333 positive**; v0.5 other-real-plant comparison **147/333 positive**. Transitions: 132 positive→negative; 10 negative→positive; 137 stayed positive; 54 stayed negative.
- Group-level directions are not driven by one arbitrarily poor alternative set. Each of the eight sets has actual-host positive in **145–165/333** species; their paired median fractions are slightly negative. Across 333 species, 129 have at least five of eight plant groups favorable to actual hosts; 158 have no more than three favorable.
- Median signed modern current-climate mismatch `actual - alternative` is **−0.01468**; median absolute mismatch gap **0.03122**. Modern climatic matching was **not exact**.

### Sensitivity to modern-climate imbalance (post hoc)

| Absolute gap cutoff | Number of species | Actual known hosts larger | Median paired fraction difference |
|---|---:|---:|---:|
| <=0.01 | 55 | 30 (54.5%) | +0.0815 percentage points |
| <=0.02 | 118 | 55 (46.6%) | −0.0950 percentage points |
| <=0.05 | 226 | 105 (46.5%) | −0.0997 percentage points |
| All | 333 | 147 (44.1%) | −0.1685 percentage points |

These cutoffs were inspected after the v0.5 scores and are **descriptive only**. The small 0.01 subset does not justify retuning to claim support.

A 4,000-draw *descriptive* bootstrap resampling **247 insect genera** (not reference draws) yields a median paired difference interval **[−0.2921,−0.0431] percentage points**; fraction-positive **[39.0%,49.3%]**. These are **not confirmatory uncertainty intervals**: alternative plant taxa are reused across insects and can induce dependence across genera.

## What is supported, and what remains open

**Supported at the screening-method level:**

1. The CHELSA-based host-history proxy contains measurable species-dependent variation beyond current conditions; across the full externally supported domain its median is small (3.7%).
2. The v0.4 apparent excess over same-insect pseudo-host sampling does not establish host-specific ecology.
3. Once the alternative is an actual **other plant species** with a matched modern climatic envelope, there is **no systematic positive specific-host advantage** (v0.5). The empirical direction is modestly opposite to a positive difference.

**Not supported by these data:**

- That glaciers left genetic barriers through insect-specific larval host resource corridors.
- That modern plant occurrences identify LGM occupied ranges, refugia, dispersal routes or reproductive isolation.
- That the original LGM genetic confirmatory programme has been reopened or rescued. Its **NOT_EVALUABLE_EMPIRICAL_RESPONSE** remains closed.

**Next biological discriminator:** leave this climate-proxy result frozen. For a genuinely new ecology programme, reconstruct **independently native-range screened, species-specific host plant distributions at LGM** with explicit land barriers and dispersal limitations, then derive *resource graph connectivity*, not only climate similarity along five quadrature points. Test on truly fresh, source-disjoint genetic sampling geometry and response panel with qualifications frozen before opening genotypes. If such fresh data cannot be assembled, report this as a climate-history comparative macroecology study **without claiming host-mediated genetic divergence**. Ecology implementation belongs in a focused biology repository; TTF retains source and inferential audits.

## Verified source and audit receipts

- v0.5 model climate workflow [#39233414704](https://github.com/zuizui0223/TTF/actions/runs/39233414704), **SUCCESS**.
- The frozen 335-insect/8-group 0BP-only plant control selection in [#39231934731](https://github.com/zuizui0223/TTF/actions/runs/39231934731), **SUCCESS**, archived on PR #149.
- Independent source-bound v0.4/v0.5 comparison [#39247018567](https://github.com/zuizui0223/TTF/actions/runs/39247018567), **SUCCESS**: `benchmarks/exploratory/lgm_host_identity_audit_v05/result.json`, receipt SHA bound to original CSVs.
- v0.4 original `species_same_niche_null_v0.4.csv` SHA256 `678c513b015af3514d2d54b5bbeea78f02a73d7c00e70c25276ce27ea0cce0a2`.
- v0.5 `species_v05.csv` SHA256 `591a1383b739d73c8ffc2cbae5c220dae4b7730629e6d447be80542b88124446`; `plant_group_v05.csv` SHA256 `236931f93d2efcfebb6988d87291ed6628e63a6ce0ae2cb227cc0825d637cf73`.
- This work remains a stacked **draft PR #149**, not merged. The repository-wide legacy CI errors are not silently reclassified as passing tests.
