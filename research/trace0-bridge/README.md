# TTF-Q Stage 13 — exact original 0 BP climate sampling bridge

**Date:** 2026-10-29. **Scope:** engineering to unblock a previously impossible empirical same-product current-climate comparator. This does NOT reopen relational v0.3 or inspect DNA/genetic responses.

## What is new

The 653-species frozen TTF-Q C1 ecological panel (96,607 selected occurrences, development 15,625 dyads, confirmatory 40,602 dyads) has been recovered and verified. A new GitHub Actions workflow retrieves the exact four original CHELSA-TraCE21k v1.0 0 BP raster files on GitHub-hosted runners, validates each file's original C1 SHA-256 **and byte size**, samples exactly the 96,607 frozen locations, and uploads only a small 96,607×4 climate-value array with a SHA-256 receipt. This removes the requirement to download 2.334 GB of rasters into the chat runtime. Source coordinates are transferred to GitHub as an unlabeled, LZMA/base64-encoded sequence without species names or GBIF record IDs. That encoded dataset remains publicly decodeable on a public branch; do not use this bridge for private or sensitive occurrence coordinates.

The on-disk value extractor is `github/sample_trace0_exact.py`. Its GitHub branch representation is `research/trace0-bridge/sample_trace0_exact.py`; the matching workflow is `.github/workflows/ttf-q-trace0-exact-sampling.yml`.

## Frozen fidelity preflight — completed

| Check | Observed |
|---|---|
| frozen climate-admissible panel | 653 species; 96,607 occurrences |
| original development dyads | 15,625 |
| original confirmatory dyads | 40,602 |
| coordinate byte SHA256 | `ae5d5b1e192500e7b12b42789b4e8ee7931f439e8c46905b267b0554ef25cef6` |
| frozen original C2 historical unique fraction | development `0.29893812458713254`; confirmatory `0.3163113642780213` |
| raw 0 BP climate asset sizes | 527,708,089 + 173,707,693 + 400,534,930 + 1,232,526,033 = 2,334,476,745 bytes |

No new 0 BP ecological result is claimed until the GitHub Actions artifact returns and its four input checksums agree with the original C1 receipt.

## Reproduce local preflight

```sh
python scripts/test_stage13_preflight.py
```

## Remote acquisition / sample bridge

The tracked workflow is triggered only by a push to the dedicated branch `research/ttf-q-trace0-exact-v01` affecting these bridge paths or via explicit `workflow_dispatch` when available. Job: download 4 public rasters, SHA-256 verify, sample all frozen coordinates, upload an artifact `ttf-q-trace0-exact-samples` with:

- `trace0_climate_values.npy` — four raw climate columns in exactly the frozen row order;
- `trace0_exact_samples_receipt.json` — C1 coordinate hash, 4 original raster hashes and sizes, values-matrix hash, schema.

No source species identity, genetic responses, or unprocessed rasters are uploaded in the artifact. The bridge's coordinate payload in the public branch is *unlabeled* but not confidential.

Once the completed artifact is available, place both its files in a `samples` folder and run:

```sh
python scripts/evaluate_trace0_exact_samples.py \
  --c1 inputs/c1_archive.zip --c2 inputs/c2_archive.zip \
  --samples samples/trace0_climate_values.npy \
  --receipt samples/trace0_exact_samples_receipt.json \
  --output results/same_product_actual
```

`results/same_product_actual/ttf_q_trace0_same_product_comparison.json` will provide for each of the two frozen panels:

1. original CHELSA V2.1 1981–2010 current-climate controlled historical unique fraction (frozen C2 reference);
2. CHELSA-TraCE21k 0 BP **same-product** current-climate controlled historical unique fraction;
3. the historical unique fraction after controlling BOTH current representations;
4. between-product pairwise current-climate overlap correlation and absolute difference.

Identical species/dyads/coordinates/PCA-KDE-Schoener-D relation formula are required. Each missing/corrupt raster or changed sample row aborts the calculation. No post-result predictor selection is permitted.

## Limitations and interpretation

**This is an empirical-data acquisition workflow, not a new observed biological result.** Current climatic similarity uses fixed 0 BP fields at *today's sampled occurrences*, while the historical relation uses the 21 ka minus 0 BP difference at those same coordinates. The historical relation is not an inferred LGM species distribution. A positive unique fraction is a predictor-information statistic, not evidence of genetic historical memory. Cross-product comparison controls the known CHELSA V2.1 vs TraCE21k mismatch but cannot exclude errors common to the single paleoclimate model.

When the same-product comparator runs, keep both existing TTF-Q panels and report them both. If the residual remains sizable, the next independent test must connect it to a separately held-out biological outcome rather than repeating predictor-only adjustments.
