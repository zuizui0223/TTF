# Genetic TTF development goal

## Scientific goal

Determine whether a spatial field learned from training species predicts within-species genetic differentiation in unseen species after the frozen endpoint-safe IBD adjustment. A qualified test can return significant transfer, qualified non-significant transfer with or without self support, or an explicit non-evaluable endpoint. The goal is an interpretable result, not a required positive result.

Predictive transfer does not causally partition geography from shared demographic history. A species-disjoint split is not a phylogenetically independent split.

## Completion criteria

1. **Complete result interface:** the authorized Phase-4 run emits the primary post-IBD statistic and inference, the within-species self diagnostic, and the predeclared descriptive total-transfer summary on the same frozen evaluation species.
2. **Verified interpretation isolation:** the total descriptor cannot alter the primary p-value, selected references, decision branch, species split, or qualification. Constant vectors retain the inherited zero-score convention; non-finite descriptive scores do not silently change the evaluation population.
3. **Fresh empirical completion:** an untouched authenticated phylogatR archive passes the sequential Phase-1/2/3 requirements before Phase-4 identity opening. A failed required gate closes the corresponding route without a biological absence claim.

This development increment implements and tests criteria 1–2 using constructed arrays only. Criterion 3 remains blocked by the source requirement in issue #21. No empirical genetic conclusion or new qualification PASS is produced by this increment.

## Total genetic transfer: implementation contract

Parent: `docs/supporting/genetic_phylogatr_phase4_response_rule_v0.1.json`, `secondary_descriptive_estimand`. The parent rule already declares this quantity descriptive unless separately qualified. This document fixes the previously missing computational details before fresh outcomes are opened.

- Input: the same in-memory locality-pair mean p-distances as the primary analysis; no additional source opening.
- Graphs, training/evaluation species, species weights, bandwidth, prior and segment quadrature: inherit the primary prepared design exactly.
- Turnover: within-species rank of raw genetic distances, **before biological IBD residualization**.
- Training response: apply the inherited within-species edge-length-rank orthogonalization to raw turnover. This geometric control is retained, so “total” here means before biological IBD adjustment, not an entirely geometry-unadjusted surface.
- Evaluation response: raw genetic-distance rank on all frozen evaluation species.
- Summary: species-equal mean Spearman score. Ties use average ranks. Constant response or prediction vectors score zero under the inherited TTF convention.
- Non-finite species scores: serialize as JSON null, retain their identities, and set the aggregate to null with `NOT_EVALUABLE_DESCRIPTIVE`. Never average only over finite survivors. Invalid mappings, distances or execution arguments still raise an error.
- Inference: none. No p-value, significance label, confidence interval or qualification flag may be borrowed from the primary residual reference family.
- Interpretation: the descriptor can describe a numerical difference before/after IBD adjustment. A positive descriptive value is not a significant total-transfer result; subtracting the two rank correlations is not a fraction of genetic variation explained by IBD.

The runner computes the descriptor after fixing its primary decision and adds `secondary_total_genetic_transfer` to the existing result payload. Primary and self calculations remain unchanged. Sequence identities and edge-distance arrays are not serialized.

The modified scorer and runner are already in `PHASE4_CODE_PATHS`: a future Phase-4 authorization binds their exact code hashes. Existing authorizations cannot silently execute changed code. No frozen qualification threshold, null family, graph or source-access decision is revised.

## Validation

- Independent dense Gaussian integration and pairwise-counting ranks reproduce descriptive scores, including tied distances, across different execution chunk sizes.
- Primary statistic, training strength, per-species scores, IBD residuals and input arrays remain identical before and after descriptive scoring.
- All constant vectors retain all evaluation species with zero scores.
- Malformed mappings, non-finite/negative distances and shape drift are rejected.
- Simulated non-finite score output retains the full species population and writes strict JSON nulls.
- Runner integration tests vary the total descriptor between negative, positive and unavailable values across all primary/self interpretation branches; primary decisions stay unchanged.
- A rejected opening authorization prevents both source extraction and descriptive scoring.

These are software and numerical checks on synthetic fixtures, not a fresh geometry-specific Type-I/power qualification.

## Result-to-manuscript completion

The receipt-only exporter is implemented in `scripts/export_genetic_ttf_manuscript.py`; usage and validation boundaries are in [GENETIC_TTF_MANUSCRIPT_EXPORT.md](GENETIC_TTF_MANUSCRIPT_EXPORT.md). It produces the Results passage, complete species table and hashed export manifest after checking terminal-receipt consistency. Closed receipts produce a status note only. Tests cover all interpretation branches and inconsistent inputs. This completes the output/reporting development milestone, not the fresh empirical validation.

## Next empirical milestone

The next empirical task remains receipt of the untouched authenticated archive and response-blind Phase 1. Same-purpose source hunting and tuning the failed Decker design remain closed under `benchmarks/frozen/genetic_empirical_opening_state_v0.3.json`.
