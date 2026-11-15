# Senecio 2021 genomic prediction × 2025 phylogeny and intrinsic RI: exact source crosswalk

**2026-11-15 | Original-source-backed, POST-OUTCOME DESCRIPTIVE.** No second original VCF genotype opening, no new predictor or effect fitted, no cross-species TTF claim. This source join follows the 2021 P1 outcome opening and hence must not be described as prospectively frozen hypothesis confirmation.

## Evidence lineage

- **2021 morphology and genotype:** James et al. (2021), *Evolution* 75:3115–3131, DOI https://doi.org/10.1111/evo.14387. P1 preregistration, original one-shot genotype run #43693277926 and immutable aggregate are unchanged.
- **2025 neutral phylogeny, crossing data:** James et al. (2025), *Molecular Ecology*, DOI https://doi.org/10.1111/mec.70090. Original-author public repository https://github.com/MaddieEJames/James_et_al._2025_MolecularEcology, commit `fd8e8f075a72c1420267eecfd482dbc6a312caef`. Exact source blobs: `data/Phylogeny.nex` = `a7c79d4b0cccb4b628e46fdcce45f417dc487d38`; `data/CrossesDH-ReciprocalsAveraged.txt` = `ae6147df032baf492d692dea33aa03398b794d78`. The latter has **two different columns called CrossType**, parsed by declared column *position* rather than collision-prone dictionary keys.

## Matched geography from exact population-ID pairs

| 2021 Dune–Headland pair | 2025 neutral tree both tips? | 2025 source clade | Frozen P1 AUC | 2025 reciprocal F1 seed-set RI | 2025 F1 germination RI |
|---|---|---|---:|---:|---:|
| D00–H00 | No | unknown | 0.956477 | — | — |
| D01–H01 | Yes | Eastern | 0.983046 | −0.134068 | −0.007206 |
| D02–H04 | No | unknown | 0.973985 | — | — |
| D03–H02 | Yes | Eastern | 0.007525 | — | — |
| D04–H05 | Yes | Eastern | 0.082466 | −0.165019 | +0.018301 |
| D05–H06 | No | unknown | 0.361863 | — | — |
| D12–H14 | No | unknown | 0.000000 | — | — |
| D32–H12 | Yes | Southeastern | 0.345110 | −0.292928 | −0.130592 |

The RI columns are the author's **Sobel–Chen averaged parent-normalized reproductive isolation** in the 2025 paper, not new experiment counts, and not individual-level P1 phenotypes. Negative estimates mean better F1 performance relative to the author's parental comparison; they are **not negative Fst or evidence of reverse adaptation**. A crossing measure exists for exactly **3/8** P1 locality pairs; the independent nuclear phylogeny includes both population tip labels for exactly **4/8** P1 pairs.

**Important falsification of a simplistic reading:** Three 2025 independently defined *Eastern* lineage population pairs have P1 genotype scores of **0.983, 0.008, 0.082**. The coarse binary geography/neutral-clade category by itself cannot account for an invariant genomic score orientation across these eastern pairs. This says nothing yet about within-clade neutral haplotypes, fine-scale ancestry, linkage, sequencing batch effects, causal trait loci or additional history predictors.

The original 2025 paper already established that parallel Senecio phenotypes can coexist with divergent reproductive-incompatibility processes, including strong *between-clade* Dune×Headland RI. **We do not independently rediscover that mechanism.** Our additional result is the provenance-locked retrospective overlap with a previously frozen out-of-locality genomic score and its restriction to 3 directly matched intrinsic-cross comparisons; no valid cross-layer regression can be estimated from those three pairs without severe pseudoreplication.

## Structural TTF-Q denominator

The 2025 crossing-comparison table has **18 distinct population-pair rows, from only 8 distinct parent-population IDs**. A maximum-cardinality matching of the endpoint-sharing graph contains **4** pair comparisons with no population used twice. That count is a **combinatorial upper bound for endpoint-disjoint comparisons**, NOT a statistical effective-sample-size estimate or 4 independent evolution events. All crossing populations are *S. lautus*, hence **species-heldout transfer: NOT_EVALUABLE**. The original 2025 source includes only two broad neutral clades, insufficient for a robust new lineage-heldout biological generalization.

## Provenance and machine checks

- [Frozen 2021→2025 exact key crosswalk](senecio_p1_2025_ri_neutral_crosswalk_receipt_20261115.json), independent crosswalk Actions #43694454108, linked sources fixed by SHA.
- [Crosswalk replay](senecio_p1_2025_ri_phylogeny_postoutcome_crosswalk_20261115.py), [synthetic schema/ID tests](../../tests/test_senecio_p1_2025_ri_phylogeny_crosswalk.py), workflow `.github/workflows/senecio-2025-ri-crosswalk.yml`.
- [Source-only crossing graph matching](senecio_2025_ri_endpoint_graph_gate_20261115.py) with independent synthetic tests; no outcome value is needed to compute the 18/8/4 denominator.
- No raw 2021 genotype rows, parent-level 2025 crossing replicates, plant-identifiable data, or third-party spreadsheet are copied into TTF.

## Future admissible biological test

A TTF species-holdout analysis needs a panel with multiple species/independent radiations and separately measured phenotype/environment/ancestry plus genuine heldout outcomes (e.g. Timema or Costus after exact source and join gates). For the Senecio source itself, fine-scale lineage source validation and a *separately registered* neutral-history analysis could be exploratory; the current P1 heldout labels cannot be retrospectively used to flip orientation and then reported as predictive gain. The original P1 AUC=0.464 remains canonical under its unaltered rules, while the label-informed direction-free oracle AUC=0.890 is only post-outcome descriptive.
