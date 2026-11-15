# Senecio P0: exploratory locality-pair held-out phenotype signal (2026-11-15)

**Status:** independently GitHub Actions replayed from source author trait means; **NOT** a formal TTF-Q evaluation and **NOT** a genetic or species-level transfer test.

## Frozen empirical design
- Original author data: [James et al. 2021 repository](https://github.com/OB-lab/James_et_al._2021-Evolution), commit `33222710aca75f79153e4118b409e7d77787d045`, exact Git blob `400f73bf21fefb6fcb775396f256c3447c755595`.
- Evaluation unit: **eight Dune–Headland locality pairs** (16 population means), not 481 plant individuals and not eight separate species.
- Headland − Dune, nine source standardized traits, each pair held out once.
- Predict held-out pair's trait contrast by mean of seven training-pair contrasts. Compare with zero contrast; no feature tuning, no genetic predictors.
- Rule preregistration: [2026-11-15 frozen design](SENECIO_P0_PREDECLARED_LOCALITY_HOLDOUT_20261115.md) committed before any trait response was read in this branch.

## Exact exploratory results

| Trait | Held-out skill vs zero difference | Sign agreement | Held-out SSE | Zero-difference SSE |
|---|---:|---:|---:|---:|
| `VH` | 85.1% | 8/8 | 3.3571 | 22.5627 |
| `WW` | 53.4% | 7/8 | 2.7153 | 5.8241 |
| `MSA` | 59.2% | 6/8 | 3.0368 | 7.4364 |
| `MSD` | 51.6% | 7/8 | 4.6430 | 9.5931 |
| `SBA` | 96.2% | 8/8 | 0.6304 | 16.7905 |
| `AREA` | 85.2% | 8/8 | 3.3078 | 22.3068 |
| `ELONG` | 66.1% | 8/8 | 4.9898 | 14.7322 |
| `DISS` | 17.1% | 6/8 | 8.1197 | 9.7912 |
| `CIRC` | 37.7% | 7/8 | 9.6285 | 15.4641 |

All **9 of 9** source traits have positive held-out skill under this one fixed estimator. Median trait skill is **59.2%**. Four traits (`VH`, `SBA`, `AREA`, `ELONG`) have 8/8 leave-pair-out sign agreement. This is a descriptive *within-S. lautus locality* finding consistent with the original paper's phenotypic parallelism, not a new independent replication or evidence for universal transfer.

The formula is `skill=1 - SSE[heldout prediction] / SSE[zero contrast]`. Global linear re-scaling of an original trait cannot change this ratio or the sign metric, even when the source authors standardized the full dataset before publication; however approximately 11% imputation of original architectural traits and source-based feature definitions still limit independent inference.

### Why this is scientifically interesting but not sufficient
1. **TTF-P within-species:** a repeatable ecotype-associated morphology pattern is measurable across eight locality contrasts.
2. **TTF-G still untested:** source genotype samples and phenotype means are not verified as identical individuals, and the published outlier SNP lists were selected using all pairs. Reuse of the outliers in training folds would leak held-out genetic responses.
3. **TTF-RI untested:** this dataset has no matching reproductive-isolation response.
4. **CREST/MLTR history not tested:** additional historical ancestry must come from an independent, held-out-safe neutral-marker or historical source. With only eight locality units, a high-dimensional historical adjustment would be statistically unsupported.

### Verification
- [Python replay](senecio_p0_pair_holdout_20261115.py) recomputes every quantity from the pinned original 2.2-kB author file.
- [Machine-readable receipt](senecio_p0_pair_holdout_exploratory_receipt_20261115.json) holds exact numeric values.
- [Synthetic fail-closed tests](../../tests/test_senecio_p0_locality_holdout.py) check eight pair IDs, duplicate/missing ecotype and source hash.
- [Independent replay CI](../../.github/workflows/senecio-p0-exploratory.yml) checks source Git SHA and equality between the Python output and the receipt.
- No pooled nine-trait p-value, species-heldout TTF score, genomic Fst response, fitness, or causal-adaptation claim is made.
