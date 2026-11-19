# TTF Fresh original-source transplant panel — response-blind first-shot (2026-11-19)

**Branch:** `prospective/original-gps-transplant-panel-20261119` cut from `main` at `a3a70eafaf315d66246862a3a41513a2be5553cf`. This is a **separate prospective pilot**. The stopped Bontrager 44/50 program, the PR #155 Hargreaves/ten-study source audit, the held-out 211-species genetic TTF, and all scientific endpoints remain untouched.

## Source-first choice, frozen BEFORE downloading any data

The source-qualification contract in `source_only_original_dataset_contract_v01.json` contains three *new* original studies, with no response-driven selection:

| Order | Study and original repository | Published design, **not** file-certified | Source evidence requested |
|---|---|---|---|
| 1 | Albano et al. 2025, `Trifolium repens`, [Dryad DOI 10.5061/dryad.bk3j9kdqk](https://doi.org/10.5061/dryad.bk3j9kdqk) | 96 source populations, 4 common gardens, 1,858 individuals; native Europe and introduced North America | Original `Albano_et_al_2025_EM_Experiment_Data.csv` header must contain `Population, Garden, Rep, Position, ContOrigin, Latitude, Longitude`; first-shot frozen file ID 4337649, version listed 2025-09-05 |
| 2 | Cross & Eckert 2026, `Camissoniopsis cheiranthifolia`, [Dryad DOI 10.5061/dryad.tb2rbp0h6](https://doi.org/10.5061/dryad.tb2rbp0h6) | 2005 and 2018 transplant experiments, **not yet independent-cohort certified** | `latlongs.csv` (Dryad ID 4692223) and `tplantlatlongs.csv` (4692224) are separately available on the public landing page, so geographic metadata can be audited without touching the fitness file |
| 3 | Carley et al. 2025, `Clarkia` demographic adaptation experiment, [Dryad DOI 10.5061/dryad.f1vhhmh24](https://doi.org/10.5061/dryad.f1vhhmh24) | 10 documented populations, 2 transplant experimental cohorts | Zip archive includes original population IDs and UTM NAD27 coordinates; file itself **not downloaded** in this first shot |

**Crucial:** the published descriptions are *not authenticated original source-GPS observations*. Albano describes European vs North American populations in separate gardens and experimental seasons; range/season or ancestry confounding must be handled before any adaptation inference. Cross & Eckert's variable `mDD5_hist` represents observational **historical climate**, not a 21ka LGM reconstruction. Neither can be mislabeled an original LGM effect.

## Executed first-shot official transport checks

GitHub Actions [run #44538463375](https://github.com/zuizui0223/TTF/actions/runs/44538463375) executed `dryad_original_source_access_preflight_v01.py` using the pinned public Dryad file IDs. It first ran **12 synthetic/fail-closed unit tests: PASS**, then attempted both the official Dryad API file-download and public website file-stream for each geographic/source file. No biological response rows were parsed or tested.

| Original file | Official API download | Official public file stream | Genuine CSV acquired |
|---|---|---|---|
| Albano source+garden original CSV | HTTP **401** | HTTP **403** | **NO** |
| Cross & Eckert `latlongs.csv` | HTTP **401** | HTTP **403** | **NO** |
| Cross & Eckert `tplantlatlongs.csv` | HTTP **401** | HTTP **403** | **NO** |

The workflow itself reports **SUCCESS** because its correct fail-closed result is `ORIGINAL_BYTES_OR_SOURCE_SCHEMA_NOT_CONFIRMED`, not because downloads succeeded. The machine-readable JSON is uploaded as Actions artifact `ttf-fresh-dryad-original-source-only-firstshot-v01` (90-day retention) and the access attempts are also printed in the job log.

**Current actual qualification**: original GPS authenticated from bytes **0**; original cohort identity independently authenticated **0**; historical 21ka climate contrast qualified **0**; fitness/genotype and LGM response modeled **0**.

Do **not** interpret 401/403 as data missing, private, retracted or impossible to obtain. They establish only that these specific CI-runner transport paths did not deliver the original bytes. Publisher/data-repository metadata remain public. Albano's associated Zenodo DOI `10.5281/zenodo.17056843` contains **two R code files**, not a second authenticated copy of the raw geographic/fitness CSV; it is not silently substituted as data.

## Scientific admission gate / next source recovery

A future source-recovery attempt must first verify the official repository version/attachment identity (DOI, file name, file ID, size/hash when available); then verify original `population_id → latitude/longitude`, garden ID and sample/cohort IDs, no garden-as-origin proxy, uniqueness and precision at the original source-grain level; then verify that 21ka climate supplies non-collinear information beyond present climate, geography and continent before opening effects.

If authorized repository web access, institutional download or the authors supply verifiable genuine data bytes, re-run this same **source-only** audit on that exact original rather than inventing source rows from abstracts. A separate prospective hypothesis, held-out design, power calculation and family-wise correction must be frozen **before any fitness, genetics or adaptation response fitting**.

**Status:** `PUBLIC_ORIGINAL_DATASETS_IDENTIFIED__ORIGINAL_BYTES_ACCESS_BLOCKED__NO_LGM_PANEL_CERTIFIED`. Existing negative/predictor-only historical climate conclusions remain unchanged; no PR #155 science files or main changes.
