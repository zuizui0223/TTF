# Palearctic non-Lepidoptera LGM extension v0.3

## Status

**Closed at the response-blind TTF-Q gate.** The authoritative qualification run (GitHub Actions run `36585250715`) completed the frozen external-data, climate, LGM-relation and TTF-Q sequence and returned `NOT_EVALUABLE_NONLEPIDOPTERA_PALEARCTIC_LGM_RELATION`.

The subgroup genetic response was therefore **never opened**: no v0.3 pairwise `T_st`, no `beta_LGM`, and no subpanel nucleotide identity were computed. The immutable result receipt is `benchmarks/frozen/genetic_palearctic_lgm_v03_qualification_result_receipt_v0.1.json`.

The parent 211-species aggregate genetic result is historically known, so this remains a prespecified same-archive secondary ecological analysis rather than an independent confirmation.

## Biological question

> Do Palearctic terrestrial holometabolous insects outside Lepidoptera that shared more similar Last Glacial Maximum climatic opportunity show post-IBD genetic differentiation in more similar places?

The ecological predictor is **shared climatic-refugial opportunity**, not proof that two species occupied the same realized refugium.

## Why v0.3 exists

The first Palearctic census retained 44 terrestrial insects, but 26 were Lepidoptera. That would substantially overlap the independent butterfly programme, so the route was narrowed to non-Lepidoptera before an LGM relation or subgroup genetic response was computed.

v0.2 then assigned every species to both source and target roles. Before its GBIF result completed, that was replaced because cross-role shared actors are not fully represented by the intended two-way source/target cluster procedure.

v0.3 therefore freezes a deterministic **disjoint role split**.

## Frozen panel

The exact response-blind panel contains 23 species:

- 14 Coleoptera;
- 8 Hymenoptera;
- 1 Diptera.

The role hash in `docs/supporting/genetic_palearctic_lgm_subpanel_rule_v0.3.json` assigns:

- 11 source species;
- 12 target species;
- no source-target overlap;
- 132 directed source-to-target dyads before external-data missingness.

No backfill or result-dependent resplitting is allowed.

External-data continuation requires at least 20 species with at least 9 source and 9 target clusters.

## Frozen occurrence route

The exact 23-species occurrence surface is produced under the v0.2 acquisition science because occurrence selection is independent of role:

1. exact GBIF species match;
2. 2010–2026 coordinate records;
3. RESOLVE 2017 Palearctic polygon filter;
4. exact-coordinate deduplication;
5. deterministic occurrence priority;
6. 10-km thinning;
7. maximum 200 and minimum 30 retained occurrences per species.

The pre-result producer binding is
`benchmarks/frozen/genetic_palearctic_lgm_v03_occurrence_producer_binding_v0.1.json`.

The bound producer is workflow run `36529160168` at commit
`80a0cf7ac4cb672b9525ce15f7f9f5762e05af6a`.

v0.3 source/target roles are reapplied only after the bound aggregate is complete. An unresolved `REQUEST_ERROR` is technical incompleteness, not ecological missingness.

## Frozen climate assets and relation

The authoritative climate contract is
`docs/supporting/genetic_palearctic_lgm_climate_input_rule_v0.4.json`, frozen before `R_LGM` after a byte-level audit of the exact TraCE21k bio01 raster semantics.

The exact already-staged CHELSA bytes are frozen in
`benchmarks/frozen/genetic_palearctic_lgm_v03_climate_asset_binding_v0.1.json` before the bound v0.3 GBIF result, `R_LGM`, or TTF-Q result is opened.

The predictor uses:

- CHELSA V2.1 1981–2010 bio1, bio7, bio12, bio15;
- CHELSA-TraCE21k V1.0 at 21 ka BP (index -190);
- TraCE21k 0 BP (index 20) as a diagnostic comparison to CHELSA current;
- a fixed 0.25-degree present-day RESOLVE Palearctic mask;
- pooled current-climate standardization and four-axis whitening;
- species-specific Gaussian KDE climatic envelopes with Scott factor (n^{-1/8});
- projection of the same current niche envelope onto current and 21-ka climate;
- `R_LGM` = Schoener D between normalized LGM suitability surfaces.

`R_present` is mandatory nuisance information. Additional declared controls are directed geographic co-opportunity, same order, same family, and occurrence-count imbalance.

Thresholded refugium polygons may not replace the continuous primary relation.

## TTF-Q opening gate

Before any subgroup genetic response, the v0.3 relation must pass every frozen gate:

- total unique `R_LGM` variance fraction after source/target identity and controls >= 0.15;
- effective source signal breadth >= 50% of realized source clusters;
- effective target signal breadth >= 50% of realized target clusters;
- maximum single-source signal share <= 0.15;
- maximum single-target signal share <= 0.15;
- null calibration passes at private-heterogeneity amplitudes A = 0, 1, 2;
- evaluable grid MDE at A = 2 <= 0.10.

A = 3 is a reported stress diagnostic only.

The response-blind v0.4 execution amendment was frozen after the ecological relation had been characterized but before the genetic response or TTF-Q outcome was opened. It changes no threshold. It only evaluates the three deterministic necessary gates first (unique information, source signal breadth, target signal breadth) and stops immediately if any fails; synthetic calibration is run only if all three pass.

Failure closes this route as
`NOT_EVALUABLE_NONLEPIDOPTERA_PALEARCTIC_LGM_RELATION`.

## Qualification result

The frozen external-data gate retained **21 of 23 species**, giving **10 source species, 11 target species, and 110 directed dyads**. The two occurrence failures were *Alloxysta fulviceps* (25 retained records) and *Euplectrus geometricida* (28), both below the prospectively frozen minimum of 30. There were no unresolved request errors.

The response-blind LGM relation itself was computable, but it was strongly redundant with present climate: `R_LGM` and `R_present` had Pearson correlation **0.861**. After source/target fixed effects and all declared controls, only **0.13098** of the original `R_LGM` variance remained uniquely identifiable, below the frozen **0.15** gate.

Signal breadth also failed because the identifiable residual relation was too concentrated. Effective source and target breadth were numerically broad enough (7.60/10 and 7.68/11), but the largest single source carried **21.29%** and the largest single target **19.35%** of residual signal, both above the frozen **15%** ceiling.

Under the authoritative v0.4 short-circuit rule, those three deterministic failures are already sufficient to close the route. The formal outcome is therefore **NOT_EVALUABLE** without opening the genetic response and without needing synthetic calibration for the decision.

An already-running full v0.3 qualification subsequently completed its synthetic grid after v0.4 had been frozen. That non-authoritative diagnostic also failed calibration: at A=0,1,2, null rejection rates were **6.98%, 7.50%, and 7.39%**, with Wilson 95% upper bounds **7.50%, 8.03%, and 7.92%**, all above the allowed **5%** ceiling. It corroborates the stop but does not define or reopen the authoritative v0.4 decision.

Therefore the route is **not a biological null test of the LGM hypothesis**. It is a design-level conclusion: with this response-blind Palearctic non-Lepidoptera panel, LGM climatic similarity is not sufficiently independent and sufficiently distributed across source and target endpoints to justify opening the genetic response.

The authority correction that fixes this chronology is
`benchmarks/frozen/genetic_palearctic_lgm_v03_qualification_authority_correction_v0.1.json`.

## Authorized response only after PASS

A TTF-Q PASS does **not** itself open genetics. It only permits a separate one-shot empirical authorization to be frozen.

The sole authorized model is:

`T_st = source_FE + target_FE + beta_LGM * z_R_LGM + declared_controls + error`

with one-sided alternative `beta_LGM > 0`.

No role, membership, climate variable, time slice, grid, bandwidth, control, qualification threshold, or substitute historical predictor may be tuned after qualification.

## Provenance correction

The v0.2 and v0.3 census receipts contain a 59-character truncation of the exact phylogatR archive SHA in a descriptive provenance field. The authoritative rule files already contain the correct 64-character SHA, and the stored v0.3 role partition exactly recomputes from it.

The immutable old receipts are not rewritten. The correction is recorded in
`benchmarks/frozen/genetic_palearctic_lgm_archive_hash_provenance_correction_v0.1.json`.

This correction changes no membership, role, ecological relation, qualification, or genetic response.

## Executable path

`.github/workflows/palearctic-lgm-climate-ttfq.yml` is the response-blind execution path:

bound GBIF aggregate -> v0.3 role binding -> exact frozen climate assets -> climate validity gate -> `R_LGM` relation -> TTF-Q.

The workflow contains no genetic-response execution step. The authoritative qualification failed before genetics, so the route is permanently closed without computing subpanel `T_st` or `beta_LGM`.
