# Albano 2025：WorldClim 1.4 matched present vs MPI ~22ka — first actual response-blind climate diagnostic

**2026-11-20 | Draft PR #161 | No plant fitness/genetic/trait response opened | Source-only**

## Verified primary data and paleoclimate bytes

The user's original Albano 2025 Dryad CSV `b3796db37145dd32ab148f5c2105ca0646923feab153d8857d4f90a5539f080c` has **96 populations** with internally consistent reported source locations. Two populations share a WorldClim **2.5 arcmin** climate raster cell; the geographic raster analysis therefore uses **95 unique, unlabeled source cells**, not 96 independent grid values. Their precise original lat/lon, population IDs and individual fitness/genotypes are not committed.

Official [WorldClim 1.4 paleoclimate](https://worldclim.org/data/v1.4/paleo1.4.html) explicitly states that the LGM simulations are **about 22 ka**, downscaled using **WorldClim 1.4** as the bias-correction baseline. This is not exactly 21,000 BP and does not share the **1970–2000 WorldClim 2.1** baseline used by the original Albano climate PCA. For correct direct subtraction the auxiliary history diagnostic pairs a **1960–1990 WorldClim 1.4 modern** raster with the corresponding WorldClim 1.4 LGM GCM at the **same 2.5 arcmin cells**.

| Source | Official original members | Original raster verification | 95 source-cell extraction |
|---|---|---|---|
| Modern 1960–1990, WC1.4 | `bio1.bil/.hdr`, `bio12.bil/.hdr` | ZIP CRC32, sizes, local-header identity, BIL WGS84 grid and SHA256 | **95/95** |
| MPI-ESM-P ~22 ka | `melgmbi1.tif`, `melgmbi12.tif` | ZIP CRC32, sizes, source-index offsets, GeoTIFF EPSG:4326 and SHA256 | **95/95** |
| CCSM4 ~22 ka (**frozen primary**) | Archive existence/header shown; original central index/raw raster still pending | **Not yet extracted** | 0 |
| MIROC-ESM ~22 ka | Archive existence/header shown; original central index/raw raster still pending | **Not yet extracted** | 0 |

Original archive source for LGM: [WorldClim v1.4 paleoclimate](https://worldclim.org/data/v1.4/paleo1.4.html); WGS84 units are temperature **°C × 10** and precipitation **millimeters** ([WorldClim v1 formats](https://worldclim.org/data/v1.4/formats.html)).

Successful actual source-only extraction runs:
- [MPI LGM 95 cells #44752767848](https://github.com/zuizui0223/TTF/actions/runs/44752767848) — synthetic 10 tests PASS, GeoTIFF verified, 95/95 valid.
- [WC1.4 matched baseline 95 cells #44754118631](https://github.com/zuizui0223/TTF/actions/runs/44754118631) — eight synthetic tests PASS, matched baseline BIL verified, 95/95 valid.
- Both original source-only cell receipts reside in `research/fresh_transplant/source_registry`.

## Actual first-shot source-only climate signal

Predefined auxiliary metric: `delta_BIOx = LGM_BIOx - modern_BIOx`, separately for BIO1 and BIO12, at the **95 unique climate cells**. Linear control columns: **intercept + continent + source latitude + source longitude + matched modern BIO1 + matched modern BIO12**. Standardize four continuous predictors within each training fold, fit OLS, and predict the held-out historical climate change. Spatial CV is *latitude quintiles within each continent*: 5 NAM folds and 5 EUR folds, each with 9–10 independent grid cells. All 10 training designs are full rank. The climatic change is the **auxiliary prediction target**, not a measured fitness or genetic response.

| Source-only metric | BIO1 (°C) | BIO12 (mm) |
|---|---:|---:|
| Mean LGM minus 1960–1990 baseline | **−14.1095** | **−128.3474** |
| SD of original climate change | 9.2192 | 188.9214 |
| Spatial 10-fold OOF R² predicting historical change from current climate/geography | **0.8514** | **0.6249** |
| Held-out RMSE | 3.5352 | 115.0919 |
| Held-out residual SD | **3.4740** | **115.5206** |
| OOF residual variance / original climate-change variance | **0.14199 (14.2%)** | **0.37390 (37.4%)** |

Continent-specific mean historical changes: NAM −19.5489°C / −167.21 mm (n=47 distinct raster cells); EUR −8.7833°C / −90.29 mm (n=48).

[Response-blind climate diagnostic #44754816913](https://github.com/zuizui0223/TTF/actions/runs/44754816913) passed nine synthetic tests and the full source-only actual calculation. The result was **independently recomputed** in JavaScript with direct retrieval of both original climate receipts, matching Python OOF R² and residual SD to numerical precision. Saved result: `source_registry/albano2025_MPI_historical_climate_novelty_source_only_v01.json`.

**Supported positive observation:** at these cells, **a linear model using current climate, continent and geography does not completely reproduce the historical climate-change field**. This is *predictor information* in one GCM, with substantial residual structure. It is neither evidence of a unique historical effect under arbitrary nonlinear spatial models nor evidence for local adaptation, plasticity, speciation or genetic separation. OOF residuals can also reflect extrapolation and model misspecification; these diagnostics are not biological-effect estimates.

## Gates still closed

- Precommitted **primary** paleoclimate GCM was **CCSM4**, not MPI-ESM-P. The current positive MPI result is **exploratory/sensitivity-only** and cannot be silently promoted to confirmatory primary.
- The future independent historical information result must survive **CCSM4** and the other sensitivity GCM (MIROC), plus spatial/continent/model-form sensitivity with fixed inference and limited false positives.
- The original experimental plant fitness, genotype, phenotype and outcome columns remain **unopened for this program**, and the 1858 plant rows are not 1858 independent historical climate observations.
- WorldClim 1.4 modern vs LGM is matched but differs from the Albano 2025 original modern WorldClim 2.1 inputs. The 470-row Montpellier six-column signed-modern-climate offset persists in the raw data; `WeightedCD` is independently authentic and unaffected. No in-place correction of old datasets or post hoc fitness-model refit.
- The older Bontrager **44/50 STOP**, 211-species transfer negatives and all previously frozen TTF results remain unchanged. PR #161 stays **DRAFT**; `main` unchanged.
