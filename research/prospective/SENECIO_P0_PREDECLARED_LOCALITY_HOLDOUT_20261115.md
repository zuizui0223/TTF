# Senecio lautus P0: preregistered *exploratory* locality-pair phenotype transfer (2026-11-15)

**Registration before reading phenotype values in this branch.** This is **not** a confirmatory TTF-Q or genetic TTF test; it does not reopen any frozen outcome, and does not claim new cross-species transfer. The purpose is a transparent feasibility pilot for locality-pair prediction only.

## Frozen input and unit

Source: James et al. 2021 Evolution author repository `OB-lab/James_et_al._2021-Evolution`, git commit `33222710aca75f79153e4118b409e7d77787d045`, path `phenotype/data_files/traitMeansPairs.txt`, blob `400f73bf21fefb6fcb775396f256c3447c755595`.

A **pair** is one Dune and one Headland population, matched by exact `Pair` string. Expected eight pairs, two ecotypes each, 16 source population mean rows. Each pair is one evaluation unit, even though individual phenotype file contains 481 plants. No test value is considered independent simply because multiple trait columns or individuals share a pair.

Nine original standardized (already logged/scaled by source authors) morphology columns, all predeclared: `VH`, `WW`, `MSA`, `MSD`, `SBA`, `AREA`, `ELONG`, `DISS`, `CIRC`. No trait, site, or direction may be selected based on results.

## Exact frozen analysis

1. For each pair p and trait j compute `delta[p,j] = Headland[p,j] - Dune[p,j]`. No further standardization, clustering or direction reversal.
2. For eight separate leave-one-pair-out folds, fit only a **per-trait training mean difference** `mu_minus_p[j] = mean(delta[q,j] for q != p)`. Predict the heldout `delta[p,j]` by this value.
3. Benchmark is a fixed zero-difference prediction `delta_hat=0` for every heldout pair and trait. For each trait, compute `skill=1-SSE(training_mean_prediction)/SSE(zero_prediction)`. Positive skill means improvement over predicting no difference; negative means worse. If baseline denominator zero, score is undefined. Report actual held-out SSE and benchmark SSE.
4. Report each trait's heldout sign-direction agreement (number of pairs for which signs of prediction and observed delta match, excluding exact zero ties as abstentions), and all traits without outcome-driven filtering.
5. Summarize the number of traits with positive skill and the median trait skill **descriptively only**. Do not compute a pooled p-value treating nine traits or eight overlapping training folds as independent experiments.
6. Predeclare failure modes: exact Git blob mismatch, missing/duplicate pair × ecotype, missing or non-finite trait value, unexpected pair count, or malformed transcript. Fail closed, without re-tuning.
7. No genomic distances, GWAS, SNP outlier selection, evolutionary history, local adaptation or sexual isolation endpoints are used. This P0 result does NOT test biological causal adaptation or any CCOC/MLTR theorem.

## Claim ceiling and next gate

A source paper showed phenotypic parallelism; P0 is an **independent re-scoring of the same published measurements**, not independent field replication. The author reports approximately 11% individual trait imputation, and published population mean traits are transformed. A locality-pair transfer summary may suggest which trait responses are repeatable but cannot certify transfer to an unseen species or genome module.

For TTF cross-lineage application, require multiple species or independently evolved lineages with non-overlapping train/evaluation units, and an independently reconstructed history/genomic predictor. Keep any new analysis in this development PR; the frozen TTF-Q main manuscript is untouched.
