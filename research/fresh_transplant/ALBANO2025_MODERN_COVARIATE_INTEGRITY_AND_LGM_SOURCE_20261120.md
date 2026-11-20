# Albano 2025：真正原本の現代気候差分不整合・主指標保全・LGM前処理資格

**2026-11-20 / 新規Draft PR #161 / response-blind／旧Bontrager 44/50停止は維持**

## 1. Original data and independent verification

User-supplied original `doi_10_5061_dryad_bk3j9kdqk__v20250905.zip`, direct content:
- `Albano_et_al_2025_EM_Experiment_Data.csv`: SHA-256 `b3796db37145dd32ab148f5c2105ca0646923feab153d8857d4f90a5539f080c` = official Dryad file ID 4337649.
- `Albano_et_al_2025_EM_WorldClim_z-standardized.csv`: SHA-256 `5adec5e8fb84b60f8070be7e20a8036e087694a1a5d1341947dcc40266945b0d`. Original WorldClim v2.1 1970–2000, six z-standardized BIO1,5,6,12,13,14 columns; 96 source rows then four garden rows.
- `Albano_et_al_2025_EM_Climatic_Distance_PCA_Scores.csv`: SHA-256 `d5866d9f217f762b9f22c8d19b8ce94ab95dd07cbf06bd7c799e6c002f58d132`. Six PCA dimensions for the same 100 sites; first three account for 97.00098% of variance.

The source file README specifies that 47 NAM rows precede 49 EUR rows, followed by Lafayette, Mississauga, Montpellier, Uppsala. **94/96 source-row correspondences are independently corroborated through a non-Montpellier common-garden BIO difference; 2 populations only appear at Montpellier and rely on README row order plus source-climate equality.** All original source IDs are retained; no biological outcome columns were analyzed.

## 2. Unequivocal Montpellier data-coordinate discrepancy (six signed source-minus-garden BIO columns)

For each original individual row and BIO variable `v` we independently reconstructed `correct_v = original_source_worldclim_z(v) - respective_garden_worldclim_z(v)`. Then compare the original `bioVdiff` column.

| Garden | Original plant rows | All 6 source–garden difference columns match? |
|---|---:|---|
| Lafayette | 465 | Yes; residual maximum < 10^-9 |
| Mississauga | 465 | Yes; residual maximum < 10^-9 |
| **Montpellier** | **470** | **No; 470 × 6 = 2,820 entries affected** |
| Uppsala | 458 | Yes; residual maximum < 10^-9 |

For all Montpellier plants and all six variables, **the published `bioVdiff` equals the source's standardized climate value, with the garden subtraction omitted**. The error `published bioVdiff - correct difference` is the *same constant per variable across all 470 rows*:

| BIO | Constant Montpellier offset (z units) |
|---|---:|
| BIO1 | +0.811598694 |
| BIO5 | +0.505914571 |
| BIO6 | +0.917345957 |
| BIO12 | −0.954782152 |
| BIO13 | −0.306477819 |
| BIO14 | −1.447067214 |

No coordinate, source ID, garden ID, PCA score, or original raw file was altered. A corrected *local-only covariate* CSV with `Garden, Position, Population` plus six reconstructed fields has been created in the user runtime, and is not added to the public Git repository.

This is a **source-only covariate consistency finding, not evidence of an erroneous adaptation conclusion**.

## 3. Most important scientific safeguard: published WeightedCD is intact

The official Zenodo [original analysis R script](https://zenodo.org/records/17056843/files/Albano_et_al_2025_EM_Analysis.R?download=1) shows that **Q1–Q3 primary analyses use `WeightedCD`**, while the Q4 per-garden six-BIO models use individual `bioVdiff + I(bioVdiff^2)` covariates. All of these are coded predictors; their relation to biological outcome values was **not** re-estimated in TTF.

We independently calculated `WeightedCD = sqrt(sum_{k=1}^3 w_k (PC_source,k - PC_garden,k)^2)` from the original 100-row PCA table, using `w = [0.6578950731, 0.2207720341, 0.0913426988]` (explained-variance proportions over *all six components*, not renormalized over the top three). **Every one of the 1,858 plants matches the published WeightedCD within 1.01 × 10^-9**, including Montpellier.

Therefore **do not claim the paper's main climatic-distance evidence is invalidated by this finding**.

For within-Montpellier Q4 models with an implicit intercept, `x_raw=x_true+c`; a full `α+βx+γx²` model is algebraically invariant under this constant shift after reparametrizing its intercept and linear coefficient:

`α_raw+β_raw(x_true+c)+γ(x_true+c)^2 = (α_raw+β_raw c+γc²)+(β_raw+2γc)x_true+γx_true²`.

Consequently in an otherwise identical model, fitted values and the quadratic shape are invariant, but the **numeric linear coefficient, its individual significance and the interpretation of climate distance = 0** need care. No original mixed-model fit/likelihood was executed and no p value or fitness result was reinterpreted. Q4 models differ from nonlinear `abs(x)`, quadratic-only-without-linear-term, or other distinct model families.

## 4. Raster grain: 95 nominal 2.5' source-cell locations (not 96)

The genuine 96 population latitude/longitude pairs are unique at 30 arcsec. For a theoretical global grid aligned with longitude −180, latitude 90:
- 30 arcsec: **96 distinct cells**
- 2.5 arcmin: **95 distinct cells** (one pair of source populations shares a cell)
- 5 arcmin: **94**
- 10 arcmin: **93**

These are **nominal geographic bins**, not values extracted from an official raster. Exact WorldClim raster origin/nodata, coastlines, temporal climate history, GPS survey precision and independent sources of stochastic information remain unverified.

## 5. Precommitted historical climate source and verified access (NOT extraction)

To avoid incorrectly mixing historical period products, `lgm_paired_climate_source_only_contract_v01.json` specifies:
- **WorldClim v1.4 matching modern baseline 1960–1990** at 2.5 arcmin (not the original Albano v2.1 1970–2000 baseline);
- **WorldClim v1.4 downscaled CMIP5 LGM snapshot ~22,000 years BP**, not a precise 21,000-year observation;
- primary CCSM4 GCM and fixed sensitivity MIROC-ESM + MPI-ESM-P; six BIO variables above and BIO1 / BIO12 primary feature emphasis;
- source-level continent, modern climate and geographic controls, with spatially blocked validation and separate new holdout specification before any phenotypic outcome analysis.

Official source pages:
- https://worldclim.org/data/v1.4/worldclim14.html
- https://worldclim.org/data/v1.4/paleo1.4.html

[GitHub Actions #44747950778](https://github.com/zuizui0223/TTF/actions/runs/44747950778) executed **10 source-only regression tests PASS** and four GET-Range preflight requests. All **four official archives replied HTTP 206 with a valid first 512-byte ZIP signature**. Reported archive sizes: paired modern 129,319,755 B; LGM CCSM4 219,610,901 B; MIROC-ESM 219,831,681 B; MPI-ESM-P 223,064,411 B. Their 512-byte prefix hashes are frozen in `source_registry/worldclim_lgm_2_5m_archive_prefix_receipt_v01.json`.

**No full archive has been downloaded/hash-validated; no LGM climate values were extracted; no GCM residual independent signal, p value, phenotype or adaptive effect has been tested.** This matters because old transferability or LGM claims must not be rescued by a post-hoc interpretation of environmental covariates.

## 6. Research decision and next operational gate

**Now confirmed:** Original Albano CSV authentic; 96 population-specific reported coordinates, 1,858 original plants, 4 gardens, **Montpellier 470-row × 6-column covariate offset**, **WeightedCD intact**, and actual WorldClim matched source endpoints accessible on CI. Reproduction code in `albano2025_modern_climate_difference_integrity_v01.py` and regression tests; original local source-only correction output in the user session. GitHub records only aggregate receipts and source code, **never raw individual source coordinates or the original outcome CSV**.

**Still HOLD:** independent historical climate residual after modern climate/geography/continent, model disagreement at ~22ka, source GPS surveying provenance, and historical-climate-to-fitness effects. Biological result remains undetermined. PR #161 is Draft and `main` unmodified.
