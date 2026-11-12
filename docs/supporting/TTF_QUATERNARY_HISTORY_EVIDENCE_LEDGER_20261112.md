# Quaternary climate evidence ledger — 2026-11-12

This working-branch ledger is a claim-boundary summary, NOT a change to frozen
TTF-Q scientific results, one-shot genetic responses, or main. All endpoints
must remain in their originally frozen families.

## Supported positive results (ecological predictors, not genetic effects)

| Evidence | Exact empirical or source-only measurement | What the positive establishes |
|---|---|---|
| TTF-Q 653 species | 15,625 development + 40,602 confirmatory directed dyads; after present climate, taxonomy, sampling and geographic opportunity, about 82–83% of the initially unique historical predictor information survives geographic adjustment | Present-climate similarity cannot simply substitute historical displacement information in the paired covariate geometry |
| Same-CHELSA-product 0BP bridge | Unique history fraction 0.3075225 / 0.3229562 (development/confirmatory); with both present-climate comparators 0.2946036 / 0.3140511 | Historical predictor nonredundancy is not explained by mixing CHELSA V2.1 with CHELSA-TraCE21k 0BP |
| External host/insect 397-species climate comparison | Median unique historical-host predictor fraction 0.0371416 after self history, two present mismatch controls and log geographic distance; 81/397 exceed 0.10 and 37/397 exceed 0.20 | Small but nonzero descriptive unique predictor fractions persist under independent, spatially separated reference sampling |
| North American mammals 80-species Stage21 | Predictor-only training historical residual unique fraction 0.6512147 after 11 controls; 5/200 false positives under synthetic spatial null, 200/200 detections for artificial positive log-odds 1.0 | Conditional known-truth gate passed on fixed inputs, NOT a real genetic-history result |

**Not directly comparable:** 0.30 / 0.037 / 0.65 use different predictor
definitions, species sets and nuisance spaces. They must NOT be pooled as
biological explanatory power or genetic R-squared.

## Robust negative controls and biological outcomes

- The early 90-insect historical fraction 0.1800108 dropped to 0.0401520
  with independent equal-count external observations and to 0.0377527
  after geographic distance. The 18% number is not a robust effect.
- A 347-insect known-host versus same-insect-niche surrogate comparison had
  median paired +1.5701 percentage points, but the modern niche matching
  distributions were not balanced.
- With 333 insects and real alternative *plant taxa* matched on present climate,
  actual-minus-alternative historical information median was -0.0016850,
  positive in only 147/333. Host-species *identity specificity* did NOT
  survive this relevant control; do not claim host-specific LGM corridors.
- Original 211-species cross-lineage post-IBD genetic transfer:
  T=0.03257, p=0.7393. Frozen species self-structure detectability
  p=0.00899, but its uncalibrated raw self correlation is negative.
- Codistributed 326-insect fresh test: beta_G=-0.004447, p=0.969.
  Sharing the same landscape did not produce transferable spatial genetics.
- Fresh independent GenDivRange He: 687 prequalified species;
  heldout prediction relative gain +0.00016164 (+0.0162%);
  species-cluster 95% interval [-0.00148644,+0.00183419].
- Fresh CaliPopGen Fst: 39 eligible species; heldout prediction
  relative gain -0.0537718 (-5.38%);
  species-cluster 95% interval [-0.1650233,+0.0068232].
- Do not post-hoc tune, reopen or combine these distinct frozen endpoints.

## Executed next gate: Jensen mammal Stage22 predictor coverage

* Fresh independent endpoint, not a revival of He, Fst or old TTF response.
* 80 species fixed as 44 source/training and 36 wholly held-out species.
  Source case sites: 10,013; fixed surveyed backgrounds: 66,687; 76,700 total.
* Stage21 geometric 50-km source-uncertainty check already passed for all
  10,013 case sites, median displacement 35.0014 km, P95 48.6970 km.
  Geometric movement alone does not establish climate-predictor sensitivity.
* New Stage22 code on PR #153 reproduces original coordinate/label SHA, applies
  the identical deterministic within-survey 50-km jitter, holds all background
  coordinates fixed, extracts the same original SHA-verified 8 CHELSA and
  5 geographic-barrier fields AT THE MOVED CASE POSITIONS, and reports only
  aggregate completeness, per-species eligibility and source/numeric hashes.
* This stage reads **no numeric genetic response** and computes **no real
  AUC or fitted genetic effect**. It cannot by itself authorize empirical
  response opening. Stop if source SHA, role counts, minimum case/background
  support, or 50-km jitter controls fail.
* Original source author fitted barrier CSVs were not published, so geographic
  controls are independently source-reconstructed, not byte-exact reproduction.
* NEVER upload IUCN sampled grids, exact cases/jitter coordinates, source
  polygons, UNEP restricted source rasters or site-level feature matrices.
* Downstream after Stage22: independently inspect positional predictor
  missingness and source-numeric repeatability; only a fully qualified frozen
  source panel can proceed to the predeclared single species-disjoint
  break-location model with 98.33% bootstrap decision threshold. Model
  outcome is still NOT EVALUATED.

## Source of record

- TTF-Q main: https://github.com/zuizui0223/TTF
- Exact same-product current climate control: https://github.com/zuizui0223/TTF/pull/150
- Host-climate external references: https://github.com/zuizui0223/TTF/pull/147
- Same-niche and matched plant controls:
  https://github.com/zuizui0223/TTF/pull/148
  https://github.com/zuizui0223/TTF/pull/149
- Genetic cross-geography null: https://github.com/zuizui0223/TTF/pull/135
- Fresh He and Fst genetic nonpositive endpoints:
  https://github.com/zuizui0223/TTF/pull/151
  https://github.com/zuizui0223/TTF/pull/152
- New 80-species mammal prospective break-location gate:
  https://github.com/zuizui0223/TTF/pull/153

## Scientific conclusion

**Historical climatic geometry is information-rich beyond current climate
and geographic closeness, but a reproducible species-general genetic effect
has not been demonstrated.** The biological target now is WHERE phylogeographic
breaks occur, with modern climate, topography, barriers and positional
uncertainty controlled and entirely unseen species evaluated.


## Actual Stage22 execution verdict — 2026-11-12

**STOP BEFORE JITTERED CLIMATE EXTRACTION; no genetic response opened.**

- Synthetic-only Stage22 determinism/containment tests and original 80-species
  panel ZIP identity passed. Species counts, roles, and original case-label SHA
  are identical to frozen Stage20.
- Strict raw-coordinate identity failed for **all 80 species** before any
  climate/geo extraction at jittered points. Exact frozen original aggregate
  coordinate SHA: fad49f6ec345ef45c1387ef76c8b3b724caf5f39322f8e9c3da56227a7a40bfd.
  New regeneration SHA: d49951247ef2cfad76f1a19406e11f5e732cc6b21a0c61c307cd85726b151dba.
- Old 2026-10-29 original and two old 2026-10-31 independent source replays
  shared the exact frozen coordinate SHA. Two 2026-11-12 current regenerations
  reproduce their new SHA identically despite Python 3.12.16 vs exact old
  Python 3.12.15 and same reported Numpy/Shapely/pyproj versions.
- Therefore **not** a Python patch-version effect. Remaining cause and actual
  location displacement unresolved. Do **not** waive SHA or call this an
  empirical biological null.
- Evidence receipt:
  research/break-location-na-mammals/stage22_coordinate_replay_stop_20261112.json
  and GitHub Actions 43128142339 / 43128609308.
- To proceed only the exact original restricted geometry can be reconstructed
  under appropriate permissions and must reproduce the frozen SHA; otherwise
  the new positional sensitivity branch ends as NOT_EVALUABLE_POSITIONAL_REPLAY.
  The older Stage21 synthetic pass and He/Fst/COI outcomes remain unchanged.


## Stage22 and Stage23 response-blind resolution (2026-11-12)

- Original October 2026 frozen coordinate SHA was obtained on Ubuntu 24.04.
  November ubuntu-latest used Ubuntu 26.04 and failed 80 of 80 coordinate
  identity checks. Explicit Ubuntu 24.04 recovered the frozen coordinate SHA
  for all 80 species. The exact low-level numeric cause is not isolated.
- Stage22 first Ubuntu24 success: 76,700 records (10,013 jittered case,
  66,687 fixed background), all 13 fields complete, all 80 species eligible
  (44 train, 36 held-out). First-run jittered 13-field numeric SHA differed
  from a later reproduction. Its raw numerical matrix was not retained.
- Two later fully independent source-only replay runs 43143255816 and
  43143333653 produced byte-identical 13-column SHA256:
  8ac29ed8c1a34dc68c245c3b1cdaeb863485d8b2f312c8881e9e03283d3a9331.
  All 13 per-column raw hashes, 3-decimal rounded hashes, numeric moments,
  14 original source assets and species support match exactly. The original
  discrepancy is preserved in the separate source-only numeric amendment.
- Stage23 run 43142984279 compared original and frozen <=50-km jittered
  CHELSA 21ka/current climate at 10,013 published case positions. Historical
  climate displacement original/moved medians 1.779670/1.773606. Pooled
  case-site Spearman 0.966726; 80-species median within-species Spearman
  0.846964, species-bootstrap descriptive 95% interval 0.786873 to 0.868038.
  Median absolute paired history shift 0.057136, P95 0.424540; 12/80 species
  rho <0.5, two negative, often in small within-species samples.
- These are ecological PREDICTOR and SOURCE POSITIONAL STABILITY results,
  NOT evidence of genetic break-history association.
- One-shot actual heldout 44/36 mammal break-location AUC is separately
  frozen in research/break-location-na-mammals/stage24_empirical_break_location_one_shot_contract_v01.json.
  Its workflow must verify ORIGINAL Stage21 13-column SHA
  51cf2953dda21ebf06f5cb78bbc0cf2c46f7a0d2bce1f7606128b27239bec2b8
  AND JITTER Stage22 amended 13-column SHA
  8ac29ed8c1a34dc68c245c3b1cdaeb863485d8b2f312c8881e9e03283d3a9331
  BEFORE opening any real genetic break-site AUC. No retuning, selecting
  subgroups, reopening older response, or exporting licensed source grids.
