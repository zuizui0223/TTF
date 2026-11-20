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


## 2026-11-19 second-wave: official Dryad metadata authenticity PASS

The public Dryad `GET /api/v2/datasets/{encoded DOI}/versions` and `GET /api/v2/versions/{id}/files` API worked from GitHub Actions even though the raw file download endpoints returned 401/403. The HAL version and file IDs were parsed from official `_links.self.href` and `_links.stash:download.href` rather than expecting nonexistent inline `id` fields. **First verified run:** [Actions #44540122635](https://github.com/zuizui0223/TTF/actions/runs/44540122635), API-only and 12 preflight tests success; expanded HAL-link test run [#44540180211](https://github.com/zuizui0223/TTF/actions/runs/44540180211), **15 tests PASS**.

| Original dataset DOI | Frozen official version ID | Exact official original file | Official file ID | Official size | Repository-declared SHA256 |
|---|---:|---|---:|---:|---|
| `10.5061/dryad.bk3j9kdqk` | **394175** | Albano Experiment_Data.csv (full original filename in JSON receipt) | 4337649 | 421805 | `b3796db37145dd32ab148f5c2105ca0646923feab153d8857d4f90a5539f080c` |
| `10.5061/dryad.tb2rbp0h6` | **435125** | `latlongs.csv` | 4692223 | 375 | `2f874477d3d0b50db9817af47c92b2b3b4d46a7ab60e2a6f86f55b0271d09c01` |
| `10.5061/dryad.tb2rbp0h6` | **435125** | `tplantlatlongs.csv` | 4692224 | 157 | `bf84ee8116386da56415c095c743a9320a28559ac8807912d84c49a7da5a3892` |
| `10.5061/dryad.f1vhhmh24` | **337849** | `Clarkia-LTREB-transplant-archive-R2.zip` | 3766673 | 7139083 | `de625d972c5f7bb2ad01709c1b0660099521d986b60b0a344a779c7ec5fcddaa` |

**These SHA256 values are given by the ORIGINAL REPOSITORY'S API and have NOT been recomputed by hashing downloaded file bytes.** They certify the published file manifest identity, **not** the original GPS coordinates or biological outcomes. The exact identity table is permanently stored in `source_registry/fresh_dryad_original_file_identity_v01.json`. A new fail-closed verifier `verify_frozen_official_dryad_identity_v01.py` and twelve known-truth tests reject DOI/version/file-size/SHA/ID drift on every subsequent official metadata run.

**Result:** four pinned original files in three archived studies have authoritative immutable identities; actual contents still inaccessible via tried authenticated/public file download endpoints. This is a substantive source-provenance advance, but NOT a new LGM/adaptation confirmation. New independent original source-GPS certified **0**, new original independently verified experiment cohorts **0**, 21ka historical contrasts tested **0**.

**Next operational dependency:** obtain original dataset bytes via an authorized download supported by Dryad or a valid repository/author mirror, verify the downloaded SHA256 against these frozen original-manifest digests, and select only source ID, coordinate, garden and cohort columns. Never infer raw values from an article map or substitute the related Zenodo *analysis code* for original CSV records.
