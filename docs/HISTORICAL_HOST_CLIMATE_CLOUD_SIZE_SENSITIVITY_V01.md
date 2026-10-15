# TTF: fixed host-climate reference-size sensitivity on the existing Costa Rica CHELSA 2D corridor pilot (2026-10-15)

**Exploratory after the first LGM/current climate-only corridor outcome was already seen. No COI identity, genetic distance, test p-value, confidence interval, or new confirmatory evidence is included.**

## Scientific question

The original native CHELSA four-variable least-cost landscape for *Bungalotis erythus* used **665 contemporary WGSRPD3 climate support vectors for the exact three HOSTS/WCVP recorded host plants**, but only **26 sampled insect localities** for the insect's own contemporary climate reference. Is the contrast between the two LGM minus present potential climate corridor costs robust to simply equalizing these sample counts?

**Equal sample size is NOT equal population sampling, geographic range, diet breadth or ecological niche coverage.** This explicitly tests numerical support-size sensitivity only, not genuine ancient host occupancy or coevolution.

## Prospective-to-this-exploratory-computation procedure

The source/method/seed/repetition count was committed before the 50-draw result as
`docs/supporting/historical_host_lcp_reference_size_rarefaction_contract_v0.1.json`. It is clearly **post-genetic-opening and post-CHELSA-LCP-exposure**, not a new genetic preregistration.

The exact original, unmodified `101×57` native CHELSA region, 26 localities and 64 original graph edges were used. For each of 50 deterministic PCG64 draws, select **26 of the original 665 host climate reference vectors without replacement within a draw**. Recalculate the **entire 2D 4-neighbor minimum resistance route** using original time slices -190 (LGM) and +20 (0BP), the original source-pinned four bioclimatic variables and their common 6,343-vector standardization. Preserve all endpoint pairs, including the one same-pixel zero route. The original insect self climate reference remains exactly 26 points; it is not artificially expanded to 665.

Before reading results, the original full host 665-point and insect self 26-point current/LGM **four 64-element arrays were reproduced with maximum absolute difference 0** from the original 2026-10-13 historical LCP CSV, SHA256 `7d255b208a485e7bc6d02e4e5748781c5e10170ed8bf161542dce3135c897b3b`.

Run: [#38675328114](https://github.com/zuizui0223/TTF/actions/runs/38675328114), **SUCCESS** (synthetic no-genotype tests and exact original CHELSA source/hash checks). Original ZIP raster SHA256 `0a797c609d26a185ba08d209a71203a584fa60f1045c4a8c335934372cee30b2`; first-layer variables `bio01,bio07,bio12,bio15`; periods -190/20. Source code `scripts/audit_host_reference_size_rarefaction.py`, tests `tests/test_hhm_cloud_rarefaction.py`.

## Results: original 64 geographical edges (climate cost-km)

| Modern climate analogue reference | Reference vectors | Mean LGM−present route cost |
|---|---:|---:|
| Frozen host climate cloud (unchanged baseline) | 665 | **+1.665783** |
| Host cloud size-matched to 26, average over 50 draws | 26 | **+5.676367** |
| Frozen insect-own sampled-climate cloud (unchanged) | 26 | **+17.794187** |

The original host minus insect contrast **−16.128403** cost-km becomes **−12.117820** cost-km when averaging the 50 matched-size host draws; its absolute size is **24.87% smaller**. This is a *numerical sensitivity*, NOT fraction of a biological effect causally attributable to reference sampling.

Across 50 host draws, mean 64-edge host LGM−present values ranged **−0.693575 to +15.850281**; descriptive 10th/90th percentiles **+1.010763 to +11.477339**; **48 of 50** drew a positive host LGM−present change. **None of the 50 host draw means exceeded the fixed insect 26-point mean +17.794187.** Across draws those values are not independent phylogeographic populations; percentiles are **not confidence intervals** and 0/50 is **not a statistical p-value**.

## Interpretation and design implication

The original larger LGM climate cost for the insect's own recorded locality climate compared with its three recorded hosts is **sensitive** to the 665:26 reference cloud point imbalance, but numerical rarefaction alone **does not reverse** the relative climate cost result in this pilot. The remaining difference could reflect real climatic niche differences *or* residual recording/geographic-sampling disparities, host breadth (three plants vs one insect), native WCVP areas as opposed to opportunistic occurrence locality positions, and model geometry.

The separate corrected NOAA ETOPO N15W090 land/illustrative static −120m assessment found **22,914/22,914 native pixels land and no static newly emergent land** within the current pilot ROI, but that cannot reconstruct ancient shorelines or ice-sheet movement. Both checks support focusing a *fresh* ecological study on historic **two-dimensional host climatic connectivity**, matched source support and genuinely time-indexed plant distributions, not post-hoc rescue of the already consumed 259-species COI hypothesis.

**The original one-shot 259-species genetic result remains `NOT_EVALUABLE_EMPIRICAL_RESPONSE`. It is not retested or substituted here.** No nucleotide identities were read in this workflow.

## Artifact identity

- Frozen results: `benchmarks/exploratory/historical_host_lcp_reference_size_26of665_result_v0.1.json`
- Exact raw source receipt SHA256: `1fd433378fb79b2fd37f3bdc6c7780eb9aede179a6144e326b1fbd415c48d435`.
- Workflow run: `38675328114`; archive artifact ID `11945442746`.
- Owner Library copy: `/TTF/historical_host_memory/host_cloud_26_size_rarefaction_20261015/`.
- No raw genotype archive, per-edge genetic responses, or genetic p-values were committed.
