# Penstemon Sierra hybrid-zone: two-zone ancestry-to-hue portability (2026-11-18)

**Status: authentic-author-source, exploratory two-zone heldout result; not cross-species TTF, not a CREST historical-memory test.** Source protocol was committed on the isolated branch before source outcome analysis. Results are exact Git source SHA verified, CI-replayed, and retained as aggregate-only JSON.

## Source

Published original research: Stone et al. (2026), *Evolution*, DOI [10.1093/evolut/qpag199](https://doi.org/10.1093/evolut/qpag199), *Penstemon davidsonii* × *P. newberryi* Sierra hybrid zones. Original author repository: [benstemon/sierra-hybrid-zone](https://github.com/benstemon/sierra-hybrid-zone), original commit `0ba4842824bcd9ddc824eedea303ad1c176fc82b`; source table `1.phenotype_curation/data/allcombined_phenotype_admixture_subpopulations_spreadsheet-sierras.csv`, original Git blob `7cd5abd0f2a5ad6de1563d3ddd2c0d857b037de8`.

Author original table: **421** distinct individuals, including **283 classified hybrids** with original continuous `modified_hue` and genome-wide `K1` ancestry coordinates both present: **Virginia Lakes 63**, **Gem Lakes 220**. Another 74 named parental plants and 64 unassigned rows were excluded by source metadata, not fitted responses.

## Fixed training/evaluation design

For each independent zone holdout: fit simple OLS `modified_hue = intercept + slope × K1` on all original classified hybrids in the **other zone**, then predict all heldout hybrids. Compare with no-`K1` benchmark that predicts only the **training** zone mean modified hue. Predictor K1 had already been jointly inferred by the authors across both zones, making this a **transductive, not fully inductive out-of-domain test**. No predicted label direction flips, externally supplied genomic outlier SNP lists, hyperparameter search, or post-hoc selected subgroup.

| Fixed fold | Training n | Heldout n | Fitted training K1 coefficient | Heldout RMSE | Training-mean baseline RMSE | Heldout SSE skill | Heldout Pearson r |
|---|---:|---:|---:|---:|---:|---:|---:|
| Virginia Lakes → Gem Lakes | 63 | 220 | +22.7444 | 4.0874 | 7.4494 | **+0.69894** | +0.8490 |
| Gem Lakes → Virginia Lakes | 220 | 63 | +21.0891 | 4.0494 | 7.9841 | **+0.74277** | +0.8722 |

Skill is `1 − heldout SSE(OLS) / heldout SSE(training-mean benchmark)`. Improvement applies to the exact published *modified_hue* scale; do not interpret values as angular degrees unless validated by the author trait transform. The two-fold equal-weighted **descriptive** mean skill is **0.72085**. There are only **two hybrid zones**, so no inferential cross-zone population standard error, p value or general transfer claim is justified.

## Valid scientific interpretation and major boundaries

**New within-source predictive result:** Genome-wide source ancestry K1 is associated with the original flower hue so consistently across both measured Sierra hybrid zones that training-zone ancestry→hue fits reduce heldout squared error relative to the training-zone mean benchmark in both directions.

**Not established:** (1) a specific anthocyanin allele causally produces hue; (2) ancestry K1 conveys an independently reconstructed historical state beyond the current genome; (3) genomic barrier architecture or reproductive isolation is predicted; (4) new **species** can be predicted; (5) any four-independent-origin generality for `Penstemon` sect. `Habroanthus`; (6) independence of the K1 feature from the heldout zone's genotypes.

The **original 2026 authors already report genetic differentiation of reproductive barriers being decoupled from floral hue divergence**, including major-hue-associated loci with different genomic-cline behaviors. This pilot must not be presented as rediscovering that biological mechanism. It is solely a source-authenticated **cross-zone predictive portability** score. The result is not compared numerically with Senecio morphology P0, genotype P1 AUC, Timema RI or Costus source eligibility because those are distinct estimands, data units and experiments.

## Provenance / replay / freeze

- [Prior frozen source-locked protocol](PENSTEMON_SIERRA_ANCESTRY_HUE_ZONE_HOLDOUT_PLAN_20261118.md), source-only pre-evaluation plan commit `84fbce458dadbf08f2d31902b298067f802f91dc`.
- [Python original author-data replay](penstemon_sierra_two_zone_ancestry_hue_20261118.py).
- [Synthetic source/holdout tests](../../tests/test_penstemon_sierra_two_zone_ancestry_hue.py).
- [Machine-readable heldout result](penstemon_sierra_ancestry_hue_two_zone_result_20261118.json).
- Independent original-source GitHub Actions [#44303269859](https://github.com/zuizui0223/TTF/actions/runs/44303269859), completed successfully. Aggregate-only artifact ID `14427096712`, SHA256 `2d61b8603df11fa37808a5c3a6564e3dcc13ad04d249ac84571e2e7297f1e39d`. Original table/genotype measurements are not redistributed.

Scientific TTF-Q registry, previous CCOC/CREST claims and main branch **remain unchanged**. This is a distinct draft research PR, not a main/scientific-unfreeze proposal.
