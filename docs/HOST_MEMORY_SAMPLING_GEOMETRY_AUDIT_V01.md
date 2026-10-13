# TTF: historical host-climate opportunity is NOT historical host corridor connectivity

**2026-10-13, explicitly exploratory, post-genetic-exposure external-predictor audit.** The original 259-species one-shot empirical test stays frozen as **`NOT_EVALUABLE_EMPIRICAL_RESPONSE`**; no new p-values, species re-selection or response opening occurred here.

## What was measured

All 641 externally complete species, 71,520 original environmental predictor edges, 17,202 original locality coordinates covering the original 642 candidate species, and the original response-blind 320/321 split were source-bound by SHA-256.

Per species, on each of its frozen graph edges, the inherited predictor `M_host = D_host_LGM − D_host_0BP` was regressed against its own climate-history predictor `M_self`, present host analog distance `D_host_0BP`, and log geographic chord distance. Independently, the *mean endpoint nearest-neighbor sampling spacing*, computed from the entire original locality point cloud, was regressed against those same controls. The ordinary squared correlation between the two residuals is a **descriptive partial R²**, not a fraction of genetic variance explained. A second diagnostic adds frozen graph degree as a sampling-topology feature.

A conventional regression degrees-of-freedom adjustment was calculated to expose small-edge graph inflation. The frozen graph edges are *not* independent because they share locality endpoints, so this is an illustrative complexity penalty only; it is NOT a valid confidence interval or p-value.

## Exact predictor-only outcome

| Equal-species summary | Development 320 | Confirmatory-role 321 | All 641 |
|---|---:|---:|---:|
| Median partial R² with endpoint nearest-neighbor spacing | 3.21% | 2.82% | **3.01%** |
| Median df-adjusted counterpart | 1.00% | 0.54% | **0.71%** |
| Median partial R² adding original graph degree | 11.48% | 10.07% | **10.85%** |
| Median df-adjusted joint counterpart | 6.62% | 5.88% | **6.17%** |

147/641 (22.9%) species have nearest-neighbor partial R² above 10%; 49/641 (7.6%) above 25%. This is **source predictor geometry dependence**, not evidence for an LGM barrier. The separate 320/321 role splits are from the previous now-outcome-exposed study, **not independent genetic validation**.

| Original graph edges per species | Species | Median raw NN partial R² | Median df-adjusted NN partial R² |
|---|---:|---:|---:|
| 14–29 | 185 | 4.72% | −1.82% |
| 30–99 | 300 | 2.92% | 0.85% |
| 100+ | 156 | 2.12% | 1.61% |

The small-graph tendency must be explicitly shown instead of calling raw partial R² a robust ecological mechanism. No p-values were generated for these descriptive statistics.

## What the five-point historical host index cannot identify

The previous `M_host` sampled nearest **modern recorded host-climate analogs** at five points on the **straight great-circle** between current insect sampling localities. The current cloud was projected against the CHELSA-TraCE21k LGM and present climate. This is a **mean climatic mismatch**. It does **not** show (i) actual larval host plants at the LGM, (ii) a connected spatial corridor, (iii) temporal persistence of host associations, (iv) realized movement, or (v) a population genetic causal barrier.


### Constructive proof that the archived five-point mean cannot identify connectivity

A **synthetic, non-biological** 9×9-grid counterexample gives two landscapes with the **same cost of 20.0** at all five straight-line interior sample positions, and thus the same five-point mean, but different off-axis routes. A four-neighbor Dijkstra scan yields a least-cost distance of **10.0** when an alternative route exists, versus **141.0** otherwise. These are arbitrary grid units, not CHELSA values or reconstructed ancient migration costs. The deterministic demonstration is `scripts/probe_host_corridor_nonidentifiability.py`, its fixture tests are `tests/test_host_corridor_counterexample.py`, and its public result is `benchmarks/exploratory/synthetic_host_corridor_nonidentifiability_v0.1.json`.

**This proves a structural non-identifiability of the existing *summary measurement*.** No significance tests, host occupancy claims or new genetic data are implied by the toy.

## Separate fresh biological question (NOT a retrospective rescue)

> When present sample spacing, intrinsic insect climate opportunity and geometric distance are accounted for, do **independently inferred LGM host-climate corridor bottlenecks** predict subsequent genetic turnover in a **new, response-unopened and taxon-disjoint genetic panel**?

The scientific distinction is **history of resource-opportunity connectivity** versus merely **observed locality sampling geometry**. It requires new measured data, rather than renaming `M_host`:

1. Derive independent host-climate suitability/analogue surfaces from CHELSA v1.0 LGM and present 4-variable rasters using exactly recorded accepted HOSTS/WCVP names. Respect time-specific land/ice masks; do not call present taxon occurrences ancient host occupancy.
2. Compute a **raster/graph-based traversable corridor or least-cost resistance** between present population coordinates, permitting alternative routes. The index must detect bottlenecks missed by the five-point straight-line average. Determine treatment of sea crossing, glaciers and cell resolution **before any future genetic response**.
3. Construct climate-only insect-self, geographic chord-distance and modern sampling-isolation controls separately. Explicitly evaluate predictor uniqueness under those controls and the effect of raster resolution/model uncertainty.
4. Freeze a **fresh species-disjoint and source-disjoint** response panel, e.g. genuinely newly obtained sequences outside the already-exposed phylogatR, 259-species, 326-species and 211-species studies. If not possible, label all analyses exploratory.
5. Establish response definiteness (segregating genetic variation/finite COI response SD) *without peeking at historical host response coefficients*, and declare an exact non-rescue terminal policy for zero-variation species. Freeze finite multiplicity and independent calibration before opening genetic outcomes.
6. Only after that, allow one registered test. Do not use the current 641-species source to derive confirmatory historical host effects or retune earlier null results.

**Feasibility audit verdict:** the climatic host predictor has a component distinguishable from simple site sampling spacing in many species, but the current archived edge means **cannot establish true historical corridor fragmentation**. The research requires reconstructing continuous, connected climatically suitable opportunity surfaces and using fresh genetic evidence.

## Literature and interpretation

- Karger et al. (2023), *Climate of the Past*, DOI: 10.5194/cp-19-439-2023: CHELSA-TraCE21k is a paleoclimate reconstruction, not host-plant occupancy.
- Bishop et al. (2025), *Molecular Ecology Resources*, DOI: 10.1111/1755-0998.14052: explicit geographic sampling design matters for landscape-genomic inference.
- van Strien, Holderegger & Van Heck (2015), *Heredity*, DOI in `https://www.nature.com/articles/hdy201462`: deme geometry/habitat topology must be distinguished from isolation by distance.
- Yannic et al. (2020), *Evolutionary Applications* DOI: 10.1111/eva.12986: modeling landscape connectivity through time needs explicit grid connectivity and dispersal assumptions; it cannot be inferred from a simple climate mismatch average.

## Reproducibility and firewall

Run `scripts/audit_host_memory_sampling_geometry.py` with **exact frozen** input files (SHA256 pinned in source). Output is `benchmarks/exploratory/host_memory_sampling_overlap_v0.1.json` and contains **only aggregate geographic/ecological predictor metrics**. Never read the COI archive, nucleotide identities, p-distances, exposed per-species coefficients or deletion outcomes in this exploratory audit.
